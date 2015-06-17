# robbyAnalyzer
Copy of the Robert Pattie UCNA analyzer

Again, tested only for ROOT 5.34, known non-functional in ROOT 5.99

06/17/2015

1) Fixed the outstanding issues with storage allocation/errors. Now runs to completion on the full 2011/2012 data set.
2) Updated plots to include titles on axes.
3) Definition of TDC Diff, relative to East and West TDC signals, have been modified in the Run.C file.  Previously, the two definitions of NTDCE were (and similarly for west)

NTDCE = (TDCE - e_tdc_cut) - TDC0    and    NTDCE = TDCW-e_tdc_cut +TDC0

where the comments indicate that the TDC0 offset it to address a binning/calculation issue by recentering the signal around 150ms without smearing.  As a results, I've changed the two NTDCE definitions to both be

NTDCE = TDCE - e_tdc_cut + TDC0

06/12/2015
At this point, most functionality has been fixed - all of the functions in 'CallAnalysisTasks()' have been uncommented-out, and most things produce pdfs.

Fixes from the base version:
0) Blinding has been turned back on! -> rtime_e, rtime_w pull from the appropriate DB entries (0,1) instead of realtime (2)

1) "Load_Histograms" and "Remove_Histograms" have been added to most CallAnalysisTasks() functions.  A (currently unfixed) issue with overfilling the analyzer with histograms was fixed by deleting the contents of btr[] and bckr[], but these were later called to produce plots in the analysis.  Load_Histograms, designed to fill the contents of a single run, was added to all loops over runs requiring btr/bckr filled, then Remove_Histograms was added to the end of these loops to the scrub the contents.  This also lead to modifications of certain loops to minimize the number of open/close commands done.

2) Selectively turn on other features (such as stat. boxes in histograms) as needed

3) Modified the 'Simple Super-Ratio' script to read the background super ratio from actual background runs.

4) Modified the asymmetry vs. E function to plot from a filled histogram

5) Modified the variable list to match the branch structure of the 2011/2012 running.

Note that there are some lines which are hardcoded, which will affect running on other machines:

a) "OctetList_Dummy" has been my test run list (which itself is a modified version of OctetList_20112012), which will pull from a run list in input_files.

b) to not actually bother solving some coding issue, I hardcoded the file-out structure in Run.C as the location of my output files (it's a soft-link to the Data2 external drive on Diablo).  Change as needed.
