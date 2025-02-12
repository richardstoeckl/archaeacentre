#' Extract PDB Assembly ID from Target Name
#'
#' This function extracts the PDB assembly ID from a given target string.
#'
#' @param target A character string containing a PDB target name with assembly information.
#' @return A data frame with columns:
#'   - `target`: The original target string.
#'   - `assembly_id`: The extracted PDB assembly ID.
#'   - `auth_asym_id`: The author-specified ID.
#' @keywords: internal
#' @importFrom stringr str_extract str_c str_to_upper
#' @importFrom tibble tibble
get_assembly_id_from_target <- function(target) {
    tibble(
        target = target,
        assembly_id = str_c(
            str_to_upper(str_extract(target, "\\w{4}")), "-",
            str_extract(target, "(?<=assembly)\\d{1,2}(?=\\.cif)")
        ),
        auth_asym_id = str_c(
            str_to_upper(str_extract(target, "\\w{4}")),
            sub(".*_(.*)", "\\1", target),
            sep = "_"
        ) %>%
            sub("-\\d+$", "", .)
    )
}

#' Retrieve Polymer Entity Identifiers and Pfam Annotations from RCSB API (batched)
#'
#' This function queries the RCSB GraphQL API to retrieve polymer entity identifiers
#' (auth_asym_id and rcsb_id) along with Pfam descriptions for a given assembly ID.
#'
#' @param assembly_ids A character vector of assembly IDs.
#' @param batch_size An integer specifying the number of assembly IDs to process in one batch (default: 400).
#' @return A data frame with columns:
#'   - `entry_id`: The PDB entry ID.
#'   - `assembly_id`: The corresponding assembly ID.
#'   - `auth_asym_id`: The author-specified asymmetry ID.
#'   - `rcsb_id`: The RCSB polymer entity instance ID.
#'   - `title`: The title of the PDB entry.
#'   - `pfam_description`: The Pfam description for the entity.
#' @importFrom httr POST accept_json content stop_for_status
#' @importFrom purrr map_dfr map_chr
#' @importFrom tibble tibble
get_polymer_info <- function(assembly_ids, batch_size = 200) {
    query <- "
    query($id: [String!]!) {
        assemblies(assembly_ids: $id) {
            rcsb_assembly_info {
                entry_id
                assembly_id
            }
            entry {
                struct {
                    title
                }
            }
            polymer_entity_instances {
                polymer_entity {
                    pfams {
                        rcsb_pfam_description
                    }
                }
                rcsb_polymer_entity_instance_container_identifiers {
                    auth_asym_id
                    rcsb_id
                }
            }
        }
    }"

    baseurl <- "https://data.rcsb.org/graphql"

    process_batch <- function(batch) {
        Sys.sleep(1)

        max_retries <- 2
        attempt <- 0
        result <- tibble()

        while (attempt < max_retries) {
            Sys.sleep(1)

            resp <- httr::POST(baseurl, httr::accept_json(),
                body = list(query = query, variables = list(id = batch)), encode = "json"
            )

            httr::stop_for_status(resp, "Access to PDB server failed")

            ret <- httr::content(resp, as = "parsed")

            if ("errors" %in% names(ret) || !"data" %in% names(ret) || length(ret$data$assemblies) == 0) {
                attempt <- attempt + 1
                if (attempt < max_retries) {
                    print(paste("Retrying batch", batch, "Attempt", attempt))
                    Sys.sleep(5)
                }
            } else {
                result <- map_dfr(ret$data$assemblies, function(assembly) {
                    rcsb_assembly_info <- assembly$rcsb_assembly_info
                    entry <- assembly$entry
                    polymer_entity_instances <- assembly$polymer_entity_instances

                    map_dfr(polymer_entity_instances, function(instance) {
                        tibble(
                            entry_id = rcsb_assembly_info$entry_id,
                            assembly_id = rcsb_assembly_info$assembly_id,
                            title = entry$struct$title,
                            auth_asym_id = instance$rcsb_polymer_entity_instance_container_identifiers$auth_asym_id,
                            rcsb_id = instance$rcsb_polymer_entity_instance_container_identifiers$rcsb_id,
                            pfam_description = ifelse(is.null(instance$polymer_entity$pfams), NA_character_,
                                paste(map_chr(instance$polymer_entity$pfams, "rcsb_pfam_description"), collapse = ", ")
                            )
                        )
                    })
                })
                break
            }
        }
        return(result)
    }

    batches <- split(assembly_ids, ceiling(seq_along(assembly_ids) / batch_size))
    map_dfr(batches, process_batch)
}

#' Retrieve Pfam Annotations for PDB Targets from RCSB API
#'
#' This function retrieves Pfam domain annotations for a given set of PDB target structures.
#' The targets are expected to be results from a local Foldseek search against the PDB database.
#' The function extracts the relevant PDB assembly and chain identifiers, queries the RCSB API,
#' and returns the Pfam descriptions associated with the identified polymer entities.
#'
#' @section Background:
#' The \href{https://www.rcsb.org/}{Research Collaboratory for Structural Bioinformatics (RCSB)} provides structural and functional
#' annotations for macromolecules stored in the Protein Data Bank (PDB). \href{https://github.com/steineggerlab/foldseek}{Foldseek} is a fast search tool
#' for comparing protein structures. When searching for similar structures in the PDB, Foldseek returns a table which contains the "target" column.
#' In the case of searches against the PDB, this target column contains the PDB filename of the hit, which is not easily interpretable.
#'
#' This function automates the extraction of Pfam domain annotations for these target PDB filenames returned by Foldseek,
#' using the \href{https://data.rcsb.org/index.html#gql-api}{RCSB GraphQL API}.
#'
#' @param targets A character vector of PDB filenames with assembly and chain information as returned by Foldseek in the "target" column.
#' @param batch_size An integer specifying the number of targets to process in one batch (default: 200).
#' The API calls to the RCSB server are made in batches to avoid overloading the server.
#' @return A data frame with columns:
#'   - `target`: The original target name from Foldseek output.
#'   - `rcsb_id`: The RCSB identifier for the matched polymer entity. Note: This uses the "label_asym_id" instead of the "auth_asym_id" used in the target name.
#'   - `title`: The title of the PDB entry associated with the entity.
#'   - `pfam_description`: The description of the Pfam family associated with the entity.
#'
#' @details
#' The function operates in the following steps:
#' 1. Extracts PDB assembly and chain IDs from the target names.
#' 2. Queries the RCSB API to retrieve corresponding polymer entity identifiers.
#' 3. Filters results to match Foldseek output.
#' 4. Retrieves Pfam annotations for the identified entities.
#' 5. Merges results into a structured data frame.
#'
#' Internally, the function calls:
#' - `get_assembly_id_from_target()`: Extracts assembly ID from PDB target name.
#' - `get_polymer_info()`: Queries the RCSB API for Pfam domain information.
#'
#' @examples
#' \dontrun{
#' targets <- c("1ABC_assembly1.cif_A", "2XYZ_assembly2.cif_B")
#' pfam_results <- get_pfam_annotation_for_targets(targets, batch_size = 400)
#' head(pfam_results)
#' }
#'
#' @importFrom stringr str_c
#' @importFrom dplyr left_join join_by select mutate filter
#' @importFrom purrr map_dfr
#' @export
get_pfam_annotation_for_targets <- function(targets, batch_size = 200) {
    # Step 1: Extract assembly IDs
    assembly_ids_df <- map_dfr(targets, get_assembly_id_from_target)
    unique_assembly_ids <- unique(assembly_ids_df$assembly_id)

    # Step 2: Get polymer entity info (batched) using the new function
    polymer_info_data <- get_polymer_info(unique_assembly_ids, batch_size)
    polymer_info_data <- polymer_info_data %>%
        mutate(auth_asym_id_for_filter = str_c(entry_id, auth_asym_id, sep = "_"))

    # Step 3: Filter results to match only relevant targets
    filtered_data <- polymer_info_data %>%
        dplyr::filter(auth_asym_id_for_filter %in% unique(assembly_ids_df$auth_asym_id))

    # Step 4: Merge results
    result <- left_join(filtered_data, assembly_ids_df, join_by(auth_asym_id_for_filter == auth_asym_id)) %>%
        select(target, rcsb_id, title, pfam_description)

    return(result)
}
