# Github workflows

Github offers access to virtual machines that run on their servers. 
This directory contains the corresponding tasks, called [workflows](https://docs.github.com/en/actions/using-workflows).
Each workflow also serves as a guidline (read the subsection job:steps in the `yml` files) to perform similar actions on your computer.

Index of workflows :
- `build_Linux.yml` : 
      Build the project on a fresh Ubuntu 22.04 virtual machine.
- `build_MacOS.yml` : 
      Build the project on a fresh MacOS virtual machine, using [MacPorts](https://www.macports.org) to install MPFR.
- `build_Windows.yml` : 
      Build the project on a fresh Windows virtual machine, using [vcpkg](https://vcpkg.io/) to install MPFR.
- `documentation.yml` :
      Build the documentation with Doxygen/Latex on a fresh Ubuntu virtual machine.
      Push the PDF documentation back to github.
      Deploy the HTML + PDF documentation on the project's website.
