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
# ENZYMES.YAML
# ---------------------------------------------------------------------------------------------------------------

# Functions to parse and validate the enzymes.yaml file


def validate_enzyme_entry(entry: dict, seen_enzymes: set) -> None:
    """
    Gets an enzyme entry and checks whether its a dictionary, contains name and sequence and whether the sequence
    only consists of allowed characters, has exactly one / which is neither at the start nor end and wether the
    name is unique. Raises value errors

    Args:
        entry: dictionary, the entry containing the data of the enzyme
        seen_enzymes: set, a set of all encountered enzyme names

    Return:
        None
    """
    if not isinstance(entry, dict):
        raise ValueError("Entry must be a dictionary")

    if "enzyme_name" not in entry:
        raise ValueError("Entry is missing an enzyme name.")

    if "recognition_sequence" not in entry:
        raise ValueError("Entry is missing a recognition sequence")

    if entry.get("enzyme_name") is None:
        raise ValueError("Entry is missing an enzyme name")

    name: str = entry["enzyme_name"]
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
    Loads enzymes.yaml and returns the content as a dictionary of dictionaries, using the enzyme name as key and
    listing recognition sequence and cut index as values in a dictionary.

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
            validate_enzyme_entry(entry, seen_enzymes)
            enzyme_name: str = entry["enzyme_name"]
            recognition_sequence: str = entry["recognition_sequence"].replace("/", "")
            cut_index: int = entry["recognition_sequence"].index("/")
            enzymes_dict[enzyme_name] = {
                "recognition_sequence": recognition_sequence,
                "cut_index": cut_index,
            }

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


# ----------------------------------------------------------------------------------------------------------------
# PRIMERS.YAML
# ----------------------------------------------------------------------------------------------------------------

# functions to parse and validate the primers.yaml


def validate_primer_entry(entry: dict, seen_primers: set) -> None:
    """
    Gets a primer entry and checks whether its a dictionary, contains name, enzyme name, core sequence, selective abses
    and whether the sequences only consists of allowed characters and wether the name is unique. Raises value errors

    Args:
        entry: dictionary, the entry containing the data of the enzyme
        seen_enzymes: set, a set of all encountered enzyme names

    Return:
        None
    """
    if not isinstance(entry, dict):
        raise ValueError("Entry must be a dictionary")

    # Validating primer name

    if "primer_name" not in entry or entry.get("primer_name") is None:
        raise ValueError("Entry is missing an primer name.")

    primer_name: str = entry["primer_name"]

    if primer_name in seen_primers:
        raise ValueError(f"Entry {primer_name} is not unique.")

    seen_primers.add(primer_name)

    # Validating enzyme name

    if "enzyme_name" not in entry or entry.get("enzyme_name") is None:
        raise ValueError(f"Entry {primer_name} is missing an enzyme name")

    # Validating core sequence

    if "core_sequence" not in entry or entry.get("core_sequence") is None:
        raise ValueError(f"Entry {primer_name} is missing a core sequence")

    core_sequence: str = entry["core_sequence"]

    if not re.fullmatch(r"[ATCG]+", core_sequence):
        raise ValueError(
            f"Core sequence of {primer_name} contains invalid characters. Only A, T, C, G are allowed."
        )

    # Validating selective bases

    if "selective_bases" not in entry or entry.get("selective_bases") is None:
        raise ValueError(f"Entry {primer_name} is missing selective bases")

    selective_bases: str = entry["selective_bases"]

    if len(selective_bases) > 3:
        print(
            f"Selective bases of {primer_name} is more than three bases long. This is discouraged as selectivity degrades"
        )
        logger.warning(
            f"Selective bases of {primer_name} is more than three bases long. This is discouraged as selectivity degrades"
        )

    if not re.fullmatch(r"[ATCGNMRWSYKVHDB]+", selective_bases):
        raise ValueError(
            f"Core sequence of {primer_name} contains invalid characters. Only A, T, C, G, M, N, R, W, S, Y, K, V, H, D, B are allowed."
        )


def load_primers_yaml(filepath: str) -> dict[str, dict[str, str | int]]:
    """
    Loads primers.yaml and returns the content as a dictionary of dictionaries, using the primer name as key and
    listing restriction enzyme, core sequence and selective bases as values in a dictionary.

    Args:
        filepath: path to the enzymes.yaml file.

    Returns:
        Dictionary containing the name of the primers, the restriction enzyme, the core sequence and the
        selective bases
    """
    try:
        with open(filepath, "r") as infile:
            raw_list = yaml.safe_load(infile)

        if not isinstance(raw_list, list):
            raise ValueError(
                "YAML format error. Root element must be a list of enzyme entries"
            )

        primers_dict: dict[str, dict[str, str | int]] = {}
        seen_primers: set[str] = set()

        for entry in raw_list:
            validate_primer_entry(entry, seen_primers)
            primer_name: str = entry["primer_name"]
            enzyme_name: str = entry["enzyme_name"]
            core_sequence: str = entry["core_sequence"]
            selective_bases: str = entry["selective_bases"]
            primers_dict[primer_name] = {
                "enzyme_name": enzyme_name,
                "core_sequence": core_sequence,
                "selective_bases": selective_bases,
            }

        return primers_dict

    except FileNotFoundError:
        logger.error(f"File primers.yaml not found at {filepath}")
        raise

    except yaml.YAMLError as e:
        logger.error(f"Error parsing primers.yaml: {e}")
        raise

    except ValueError as ve:
        logger.error(f"Validation Error: {ve}")
        raise

    except Exception as e:
        logger.error(f"Unexpected error while loading primers.yaml: {e}")
        raise
