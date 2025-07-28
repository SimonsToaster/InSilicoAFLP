# ---------------------------------------------------------------------------------------------------------------
# FILE IO
# ---------------------------------------------------------------------------------------------------------------

# Functions used to read in data from files, validate their content and format it to requirements of the core
# program.

# ---------------------------------------------------------------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------------------------------------------------------------

import yaml
import re
import logging

# ---------------------------------------------------------------------------------------------------------------
# LOGGING
# ---------------------------------------------------------------------------------------------------------------

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------------------------------------------
# RESTRICTION ENZYME
# ---------------------------------------------------------------------------------------------------------------
def validate_enzyme_entry(entry: dict, seen_enzymes: set) -> None:
    """
    Gets an enzyme entry and checks whether its a dictionary, contains name and sequence and whether the sequence
    only consists of allowed characters, has exactly one / which is neither at the start nor end and wether the
    name is unique

    Args:
        entry: the entry containing the data of the enzyme
        seen_enzymes: a set of all encountered enzyme names

    Return:
        None
    """
    if not isinstance(entry, dict):
        raise ValueError("Entry must be a dictionary")

    if "name" not in entry:
        raise ValueError("Entry is missing an enzyme name.")

    if "recognition_sequence" not in entry:
        raise ValueError("Entry is missing a recognition sequence")

    if entry.get("name") is None:
        raise ValueError("Entry is missing a name")

    name: str = entry["name"]
    sequence: str = entry["recognition_sequence"]

    if entry.get("recognition_sequence") is None:
        raise ValueError(f"Recognition sequence of {name} is missing.")

    if name in seen_enzymes:
        raise ValueError(f"Entry {name} is not unique.")

    seen_enzymes.add(name)

    if sequence.count("/") == 0:
        raise ValueError(f"Recognition sequence of {name} is missing a cut site (/)")

    if sequence.count("/") > 1:
        raise ValueError(
            f"Recognition sequence of {name} contains more than one cut site (/)"
        )

    if sequence.startswith("/") or sequence.endswith("/"):
        raise ValueError(
            f"Recognition sequence of {name} starts or ends with a cut site (/)"
        )

    clean_sequence = sequence.replace("/", "")
    if not re.fullmatch(r"[ATCG]+", clean_sequence):
        raise ValueError(
            f"Recognition sequence of {name} contains invalid characters. Only A, T, C, G are allowed."
        )


def load_enzymes_yaml(filepath: str) -> dict[str, dict[str, str | int]]:
    """
    Loads enzyme.yaml and returns the content as a dictionary of dictionaries, using the name as key and listing
    recognition sequence and cut index as values in a dictionary.

    Args:
        filepath: path to the enzymes.yaml file.

    Returns:
        Dictionary containing the names fo the restriction enzyme,
        the clean recognition sequence and the cut index.
    """
    try:
        with open(filepath, "r") as infile:
            raw_list = yaml.safe_load(infile)

        if not isinstance(raw_list, list):
            raise ValueError(
                "YAML format error. Root element must be a list of enzyme entries"
            )

        enzymes_dict: dict[str, dict[str, str | int]] = {}
        seen_enzymes: set[str] = set()

        for entry in raw_list:
            print(entry)
            validate_enzyme_entry(entry, seen_enzymes)
            enzyme_name: str = entry["name"]
            recognition_sequence: str = entry["recognition_sequence"].replace("/", "")
            cut_index: int = entry["recognition_sequence"].index("/")
            enzymes_dict[enzyme_name] = {
                "recognition_sequence": recognition_sequence,
                "cut_index": cut_index,
            }
            print(enzymes_dict)

        return enzymes_dict

    except FileNotFoundError:
        logger.error(f"File enzymes.yaml not found at {filepath}")
        raise

    except yaml.YAMLError as e:
        logger.error(f"Error parsing enzymes.yaml: {e}")
        raise

    except ValueError as ve:
        logger.error(f"Validation Error: {ve}")
        raise

    except Exception as e:
        logger.error(f"Unexpected error while loading enzymes.yaml: {e}")
        raise
