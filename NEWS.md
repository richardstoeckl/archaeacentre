# archaeacentre 1.1.0

IMPROVEMENTS

* Major rework for `get_pfam_annotation_for_targets`.
    * The function is now much more resilient to unusual target filenames.
    * The function now calls the RCSB API only once leading to less strain on the server
    * The function now includes forced wait times to prevent overloading the server
    * The function now includes some fallback logic if a batch of target filenames fails on the first API call
    * The function now has the default batch size reduced to prevent overloading the server

# archaeacentre 1.0.2

IMPROVEMENTS

*  Rework batch logic for `get_pfam_annotation_for_targets` for faster processing

BUGFIXES

* Fix `get_pfam_annotation_for_targets` returning the results of other chains from the same structure

# archaeacentre 1.0.1

* Update Readme and Vignettes

# archaeacentre 1.0.0

NEW FEATURES

* Added a function to get the Pfam Annotations for Foldseek PDB target filenames

# archaeacentre 0.99.0

NEW FEATURES

* Initial release of the ArchaeaCentre package.
* Added A wrapper function around a ggplot2 that includes the default settings for microbial growth curves as used in the Archaea Centre.
