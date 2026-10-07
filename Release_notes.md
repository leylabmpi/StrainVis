## Release Notes

### Version 1.4.0

- Added executable files for Windows and Linux/MacOS.
- The requested version of Panel was changed in the strainvis.yml file to at least 1.9.3.
- Added LoadingSpinner widget when uploading a SynTracker input-file via the FileInput widget.
- Added LoadingSpinner widget to several UI elements, where the user has to wait until the plots are loaded.
- Improved responsiveness of the slider widgets.
- Updated the sample_data files (included in the release tar file).

### Version 1.4.1

- Fixed problems with the 'Show annotations' option.
- Updated the sample_data files (included in the release tar file).

### Version 1.4.2

- Added the option to execute in debug mode

### Version 1.4.3

- Added output of 'n_same_category', 'n_different_category' to the P-values tabl for export (in the 'Distribution among species' plot).
- Added a note in the initial loading page when a metadata file was uploaded.
- Added a debug version
- Launch the server dynamically using via pn.serve() instead of running the command 'panel serve...' from the command-line
  (now there is one python process for each server start instead of two).
- Removed the option to upload SynTracker input file using FileInput widget (causes problem when the file is too big).
- Fixed Windows executable files.
  