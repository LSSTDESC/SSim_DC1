# DC1 Configuration, production, validation and tools for the DC1 Data Set.

The repository holds all of the issues, tools and validations used for the DC1 production. It also holds LaTex files which document the DC1 production strategy and techniques.  There is a GitHub Issues list to track the work necessary to complete the production.  See the Issues and Milestones for more information.

A high level summary of the production plan follows:
 
- We will produce PhoSim and imSim e-image files.  
- We will also produce flats in order to correct for vignetting in the PhoSim files.
- The current CI workflow needs to be tweaked to work at NERSC.  We need to understand how to make both the workflow and DM work on full focal planes and more generally how to run in the time restricted batch environment there.
- CI and SSim will specify the cookbooks (based on Twinkles experience).  CI will run the jobs.
- The DM output catalogs will be placed in a qserv database if possible.

Last edited by C. Walter (01/26/17)
