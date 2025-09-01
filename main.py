# --------------------------------------------------------------
# IMPORTS
# --------------------------------------------------------------

import functions.file_io as io
import logging

# --------------------------------------------------------------
# LOGGING
# --------------------------------------------------------------

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    filename="logfile.log",
    filemode="w",
)

# --------------------------------------------------------------
# FILE PATHS
# --------------------------------------------------------------

path_enzymes: str = "Data\\test_enzyme.yaml"

# --------------------------------------------------------------
# MAIN
# --------------------------------------------------------------


def main() -> None:
    enzyme_dict = io.load_enzymes_yaml(path_enzymes)

    print(enzyme_dict)


if __name__ == "__main__":
    main()
