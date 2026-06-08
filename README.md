# MBL_BB_Project
Python files for data analysis and visualization for the 'Bulb Badies' team in the Neurobiology Course

## Code Installation 
To best utilize and customize this pipeline I recommend downloading and using Visual Studio Code at https://code.visualstudio.com/download.

Once Visual Studio code is installed you can clone the git repo using web URL or download a zip file containing all of the code under the 'Code' tab above. Open the repository in VSCode before proceeding.

## UV
Once the pipeline is installed on your computer you will need to download UV to ensure you have everything you need for the code to run properly.

### Installing UV
Install UV through the powershell using installation instructions outlined here: https://docs.astral.sh/uv/getting-started/installation/.

Check that UV is installed by checking the version:
`uv --version`

A common issue with installation is uv not being properly added to PATH following installation. There should be instructions that pop up, but if you get lost use instructions documented here: https://docs.astral.sh/uv/reference/installer/ 

### Installing dependencies through UV
Once UV is installed to your computer, you can install all needed dependencies with `uv sync`

After downloading the dependencies initiate a local virtual environment with `uv venv`. Try activating the virtual environment using `.\.venv\Scripts\Activate.ps1`. If you don't have permissions for this, try running in command prompt instead of powershell. To switch default terminal profile go to View->Command Pallet->Search "Terminal: Select Default Profile"->Command Prompt. Try again in command prompt.

To check that you have all necessary dependencies run `uv pip list` and make sure you have everything specified in the 'pyproject.toml' file.