"""
All service functions required to run tests
Used Hungarian name to avoid picking this module by pytest
"""

import os
import subprocess
from pathlib import Path

# === Utils for downloading test data from the s3 storage === #

TEST_DATA_FILENAME = "all_test_data.zip"
TEST_DATA_ARCHIVE_URL = f"https://wes-qc-data.cog.sanger.ac.uk/all_test_data/{TEST_DATA_FILENAME}"
TEST_DATA_PARENT_DIR_URL = "https://wes-qc-data.cog.sanger.ac.uk"
TEST_DATA_DIR_NAMES = ["control_set_small", "unit_tests", "training_sets", "resources"]


# TODO: download using a .txt file with the list of all files instead of archiving
# TODO: test this draft
def download_test_data_using_files_list(files_list: str, outdir: str) -> None:
    """`files_list` must be generated with the following command ran from the directory with the test data:
    ```
    find . ! -type d -print > ../files_list.txt
    ```
    Then it can be used as an input to this function.
    """
    with open(files_list, "r") as f:
        all_files = [file_path for file_path in f.readlines() if not file_path.startswith("#")]

    downloaded_test_dirs = [
        test_dir for test_dir in TEST_DATA_DIR_NAMES if os.path.exists(os.path.join(outdir, test_dir))
    ]
    print(f"DEBUG: Test folders {', '.join(downloaded_test_dirs)} already downloaded")
    print(f"DEBUG: Checking for file presence in {outdir}")
    for file_path in all_files:
        file_path = file_path.rstrip()
        file_path_relative_to_download_dir = file_path.replace("./", "", 1)
        data_folder = file_path_relative_to_download_dir.split("/")[0]  # TODO: optimise

        # naive approach - if parent dir of the file already exists, don't download the file
        if data_folder in downloaded_test_dirs:
            continue

        file_url = file_path.replace(".", TEST_DATA_PARENT_DIR_URL, 1)  # create download urls for each file

        file_destination_path = os.path.normpath(os.path.join(outdir, file_path_relative_to_download_dir))
        file_destination_dir = os.path.dirname(file_destination_path)
        # remove possible double slashes

        # Checking that the file exists
        if Path(file_destination_path).is_file():
            continue

        print(f"DEBUG: Downloading file: {file_url}")
        subprocess.run(
            ["wget", "-nv", "-nc", file_url, "-P", file_destination_dir]
        )  # downlaod the file into destination


def move_dirs(move_dirs: dict) -> None:
    print("Copying data to correct dirs")
    for dir_to_move, destination_dir in move_dirs.items():
        subprocess.run(["cp", "-vnr", dir_to_move, destination_dir])


# TODO: make versatile, don't download if already exists
def download_test_data_from_s3(outdir: str, move_dirs: dict, clean_up_unzip_dir: bool = False) -> None:
    """Download compressed test data from the s3 storage and
    move the subfolders to the correct destinations inside the repo.

    Parameters
    ----------
    outdir : str
        Path to download and extract compressed test data.

    move_dirs : dict
        Directories to move. Keys are the dirs to be moved, values are their destination dirs

    clean_up_unzip_dir : bool, default: False
        Remove the folder with the unzipped data after moving the data from it.

    Return
    ------
        None
    """
    print("Downloading data from the s3 bucket")  # TODO: improve logging
    # Download zipped archive with all the data
    subprocess.run(["wget", "-nc", TEST_DATA_ARCHIVE_URL, "-P", outdir])  # skip existing

    # Unzip the data
    print("Unzipping the data")
    subprocess.run(["unzip", "-n", os.path.join(outdir, TEST_DATA_FILENAME), "-d", outdir])
    # Move unzipped data to correct folders in the test dir
    print("Copying data to correct dirs")

    for dir_to_move, destination_dir in move_dirs.items():
        subprocess.run(["cp", "-vnr", dir_to_move, destination_dir])

    if clean_up_unzip_dir:
        subprocess.run(["rm", "-r", outdir])
