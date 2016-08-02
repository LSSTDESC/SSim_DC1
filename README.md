# SSim_DC1_Roadmap
Configuration, production and validation specifications for the DC1 Data Set.

This repository holds LaTex files to docuement the DC1 production
strategy and techniques.  It will likely eventually become a technical
note.  If you are working on the DC1 production please feel free to
work on the document and add yourself to the author list.

The repository also hosts the GitHub Issues list to track the work
necessary to complete the production.  See the Issues and Milestones
for more information.


A high level summary of the production plan follows:

- We will produce PhoSim e-image files.  
- We will also produce flats in order to correct for vignetting.
- The current CI workflow needs to be tweaked to work at NERSC.  We need to understand how to make both the workflow and DM work on full focal planes and more generally how to run in the time restricted batch environment there.
- CI and SSim will specify the cookbooks (based on Twinkles experience).  CI will run the jobs.
- The DM output catalogs will be placed in a qserv database if possible.

Last edited by C. Walter (8/2/16)
