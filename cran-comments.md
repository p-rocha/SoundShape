## Test environments
* local R installation, Win x86_64
* external Mac, 
* ubuntu 16.04 (on travis-ci), R 4.1.0
* win-builder (devel R 4.1.0; release R 4.0.2)

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release.

* checking for future file timestamps ... NOTE
  unable to verify current time

  This note only appeared whilst running R CMD check using R 4.0.2 on local Win x86_64. This appears to be related to R CMD check and it did not return a NOTE when I run it on R devel.
