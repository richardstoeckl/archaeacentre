#' Extract PDB Assembly ID from Target Name
#'
#' This function extracts the PDB assembly ID from a given target string.
#'
#' @param target A character string containing a PDB target name with assembly information.
#' @return A data frame with columns:
#'   - `target`: The original target string.
#'   - `assembly_id`: The extracted PDB assembly ID.
#' @keywords: internal
#' @importFrom stringr str_extract
get_assembly_id_from_target <- function(target) {
    pdb_id <- toupper(stringr::str_extract(target, "\\w{4}"))
    assembly_number <- stringr::str_extract(target, "(?<=assembly)\\d(?=\\.cif)")
    assembly_ids <- paste0(pdb_id, "-", assembly_number)
    df <- base::data.frame(target = target, assembly_id = assembly_ids)
    return(df)
}

#' Retrieve Polymer Entity Identifiers from RCSB API
#'
#' This function queries the RCSB GraphQL API to retrieve polymer entity identifiers
#' (auth_asym_id and rcsb_id) for a given assembly ID.
#'
#' @param assembly_ids A character vector of assembly IDs.
#' @return A data frame with columns:
#'   - `entry_id`: The PDB entry ID.
#'   - `assembly_id`: The corresponding assembly ID.
#'   - `auth_asym_id`: The author-specified asymmetry ID.
#'   - `rcsb_id`: The RCSB polymer entity instance ID.
#' @keywords: internal
#' @importFrom httr POST accept_json content
get_auth_asym_ids <- function(assembly_ids) {
    # Prepare GraphQL query
    query <- "
  query($id: [String!]!) {
    assemblies(assembly_ids: $id) {
      rcsb_assembly_info {
        entry_id
        assembly_id
      }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          rcsb_id
        }
      }
    }
  }"

    # Send POST request to GraphQL API
    baseurl <- "https://data.rcsb.org/graphql"
    resp <- httr::POST(
        baseurl,
        httr::accept_json(),
        body = list(query = query, variables = list(id = assembly_ids)),
        encode = "json"
    )

    # Check for errors in the response
    if (httr::http_error(resp)) {
        stop("Access to PDB server failed")
    } else {
        ret <- httr::content(resp)
    }

    if ("errors" %in% names(ret)) {
        stop("Retrieving data from PDB failed")
    }

    if (!"data" %in% names(ret) || length(ret$data$assemblies) == 0) {
        stop("No data retrieved")
    }

    # Extract and format the data
    data <- lapply(ret$data$assemblies, function(assembly) {
        entry_id <- assembly$rcsb_assembly_info$entry_id
        assembly_id <- assembly$rcsb_assembly_info$assembly_id
        lapply(assembly$polymer_entity_instances, function(instance) {
            list(
                entry_id = entry_id,
                assembly_id = assembly_id,
                auth_asym_id = instance$rcsb_polymer_entity_instance_container_identifiers$auth_asym_id,
                rcsb_id = instance$rcsb_polymer_entity_instance_container_identifiers$rcsb_id
            )
        })
    })

    # Flatten the list and create a data frame
    data <- do.call(rbind, lapply(data, function(x) do.call(rbind, x)))
    data <- data.frame(data, stringsAsFactors = FALSE)

    return(data)
}

#' Retrieve Pfam Annotations from RCSB API
#'
#' This function queries the RCSB GraphQL API to retrieve Pfam domain descriptions
#' for a given set of RCSB polymer entity instance IDs.
#'
#' @param rcsb_ids A character vector of RCSB polymer entity instance IDs.
#' @return A data frame with columns:
#'   - `rcsb_id`: The RCSB polymer entity instance ID.
#'   - `title`: The title of the PDB entry.
#'   - `pfam_description`: The Pfam description for the entity.
#' @keywords: internal
#' @importFrom httr POST accept_json content
get_pfam_annotation <- function(rcsb_ids) {
    # Prepare GraphQL query
    query <- "
  query($id: [String!]!) {
    polymer_entity_instances(instance_ids: $id) {
      rcsb_id
      polymer_entity {
        entry {
            struct {
              title
            }
        }
        pfams {
          rcsb_pfam_description
        }
      }
    }
  }"

    # Send POST request to GraphQL API
    baseurl <- "https://data.rcsb.org/graphql"
    resp <- httr::POST(
        baseurl,
        httr::accept_json(),
        body = list(query = query, variables = list(id = rcsb_ids)),
        encode = "json"
    )

    # Check for errors in the response
    if (httr::http_error(resp)) {
        stop("Access to PDB server failed")
    } else {
        ret <- httr::content(resp)
    }

    if ("errors" %in% names(ret)) {
        stop("Retrieving data from PDB failed")
    }

    if (!"data" %in% names(ret) || length(ret$data$polymer_entity_instances) == 0) {
        stop("No data retrieved")
    }

    # Extract and format the data
    data <- lapply(ret$data$polymer_entity_instances, function(instance) {
        rcsb_id <- instance$rcsb_id
        title <- instance$polymer_entity$entry$struct$title
        if (is.null(instance$polymer_entity$pfams) || length(instance$polymer_entity$pfams) == 0) {
            pfam_descriptions <- NA
        } else {
            pfam_descriptions <- sapply(instance$polymer_entity$pfams, function(pfam) {
                pfam$rcsb_pfam_description
            })
        }
        data.frame(rcsb_id = rcsb_id, title = title, pfam_description = pfam_descriptions, stringsAsFactors = FALSE)
    })

    # Combine all data frames into one
    data <- do.call(rbind, data)

    return(data)
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
#' @param batch_size An integer specifying the number of targets to process in one batch (default: 1000).
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
#' - `get_auth_asym_ids()`: Retrieves chain identifiers and RCSB polymer entity IDs.
#' - `get_pfam_annotation()`: Queries the RCSB API for Pfam domain information.
#'
#' @examples
#' \dontrun{
#' targets <- c("1ABC_assembly1.cif_A", "2XYZ_assembly2.cif_B")
#' pfam_results <- get_pfam_annotation_for_targets(targets, batch_size = 1000)
#' head(pfam_results)
#' }
#'
#' @importFrom stringr str_extract
#' @export
get_pfam_annotation_for_targets <- function(targets, batch_size = 1000) {
    # Function to process a batch of targets
    process_batch <- function(target_batch) {
        # Step 1: Use get_assembly_id_from_target() to get the assembly_id for each target
        # assembly_ids <- sapply(unique(target_batch), get_assembly_id_from_target)
        assembly_ids <- do.call(rbind, lapply(target_batch, get_assembly_id_from_target))

        # Step 2: Use get_auth_asym_ids() to get all rcsb_ids for the assemblies
        auth_asym_data <- get_auth_asym_ids(assembly_ids$assembly_id)
        auth_asym_data$auth_asym_id_for_filter <- paste(auth_asym_data$entry_id, auth_asym_data$auth_asym_id, sep = "_")

        # Step 3: Extract the auth_asym_id from the target names
        target_auth_asym_ids <- sapply(target_batch, function(target) {
            paste(toupper(stringr::str_extract(target, "\\w{4}")), sub(".*_(.*)", "\\1", target), sep = "_")
        })

        # Step 4: Filter the returned dataframe to get the rcsb_id for the auth_asym_ids that are encoded in the target names
        filtered_data <- auth_asym_data[auth_asym_data$auth_asym_id_for_filter %in% target_auth_asym_ids, ]

        # Step 5: Get the annotation for these rcsb_id using get_pfam_annotation()
        pfam_data <- get_pfam_annotation(filtered_data$rcsb_id)

        # Step 6: Merge the filtered data with the PFAM data
        result <- merge(filtered_data, pfam_data, by = "rcsb_id", all.x = TRUE, sort = FALSE)

        # Step 7: Get the rcsb_id for all targets
        result$assembly_id_for_merge <- paste(result$entry_id, result$assembly_id, sep = "-")

        # Step 7:
        result2 <- merge(result, assembly_ids, by.x = "assembly_id_for_merge", by.y = "assembly_id", all.x = TRUE, sort = FALSE)
        result2 <- result2[, c("target", "rcsb_id", "title", "pfam_description")]

        return(result2)
    }

    # Split targets into batches
    target_batches <- split(unique(targets), ceiling(seq_along(unique(targets)) / batch_size))

    # Process each batch and combine results
    results <- do.call(rbind, lapply(target_batches, process_batch))

    return(results)
}
