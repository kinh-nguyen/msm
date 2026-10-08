# Build code/msm.so once (run by code/compile.sh before the array jobs).
TMB::compile(here::here("code/msm.cpp"), flags = "-Wno-ignored-attributes", framework = "TMBad")
