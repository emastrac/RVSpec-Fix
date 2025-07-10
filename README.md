# RVSpec-Fix
Library fixes for RaveSpan to allow the program to run.
You can find the download for RaveSpan [here](https://users.camk.edu.pl/pilecki/ravespan/?page=download).

## Installation
### Windows
1. If you have a version of Python 3.9 installed, skip to step 3. Otherwise, download [Python 3.9.13](https://www.python.org/downloads/release/python-3913/) (or any version of Python 3.9).
2. Install the version of Python you just downloaded and make note of the installation location.
3. Download RaveSpan from [Dr. Pilecki's website](https://users.camk.edu.pl/pilecki/ravespan/?page=download).
4. Optionally, download the template library slib from the same page.
5. Go to the [Code](https://github.com/emastrac/RVSpec-Fix) section of this repository, and download the files from the green Code dropdown menu.
6. Extract `ravespan_p3.tgz` and `RVSpec-Fix-main.zip` (it does not matter where).
7. In the command prompt, `pip install` `matplotlib`, `numpy`, `scipy`, and `PyQt5`
* If Python 3.9 is not in your PATH, use the installation location identified to do `\path\to\Python39\Scripts\pip.exe install` instead.
8. In the command prompt, do `python3 setup.py install --user`.
* If Python 3.9 is not in your PATH, use the installation location identified 
