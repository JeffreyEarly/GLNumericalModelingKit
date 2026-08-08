# FFTW real-to-real benchmark

Status: **PASSED**

## Environment

| Field | Value |
|---|---|
| MATLAB | 26.1.0.3312084 (R2026a) Update 4 |
| Architecture | maca64 |
| Processor | Apple M5 Max |
| Memory | 48 GB |
| FFTW | fftw-3.3.8 |
| Library | /Applications/MATLAB_R2026a.app/bin/maca64/libmwfftw3.3.dylib |

## Complete-call winners

| Nz | Batch | Type | Transform | Direction | Dense (ms) | Extension (ms) | FFTW (ms) | Speedup | Winner | Eligible | Max rel. error |
|---:|---:|---|---|---|---:|---:|---:|---:|---|---|---:|
| 33 | 1 | real | cosine | forward | 0.0081 | 0.1352 | 0.1718 | 0.047x | dense-matrix | no | 3.53e-15 |
| 33 | 1 | real | cosine | inverse | 0.0084 | 0.1337 | 0.1720 | 0.049x | dense-matrix | no | 2.19e-15 |
| 33 | 1 | real | sine | forward | 0.0075 | 0.1435 | 0.1712 | 0.044x | dense-matrix | no | 6.16e-15 |
| 33 | 1 | real | sine | inverse | 0.0077 | 0.0981 | 0.1745 | 0.044x | dense-matrix | no | 4.13e-15 |
| 33 | 1 | complex | cosine | forward | 0.0075 | 0.0824 | 0.0910 | 0.083x | dense-matrix | no | 4.93e-15 |
| 33 | 1 | complex | cosine | inverse | 0.0084 | 0.1095 | 0.1020 | 0.082x | dense-matrix | no | 3.72e-15 |
| 33 | 1 | complex | sine | forward | 0.0091 | 0.1136 | 0.1809 | 0.050x | dense-matrix | no | 2.95e-15 |
| 33 | 1 | complex | sine | inverse | 0.0057 | 0.0830 | 0.1650 | 0.035x | dense-matrix | no | 4.55e-15 |
| 33 | 8320 | real | cosine | forward | 0.0912 | 1.1897 | 0.4139 | 0.220x | dense-matrix | no | 4.53e-15 |
| 33 | 8320 | real | cosine | inverse | 0.0853 | 1.3848 | 0.4167 | 0.205x | dense-matrix | no | 5.01e-15 |
| 33 | 8320 | real | sine | forward | 0.0713 | 1.4257 | 0.1975 | 0.361x | dense-matrix | no | 4.27e-15 |
| 33 | 8320 | real | sine | inverse | 0.0822 | 1.3388 | 0.4556 | 0.180x | dense-matrix | no | 4.12e-15 |
| 33 | 8320 | complex | cosine | forward | 0.3592 | 1.4250 | 0.4465 | 0.805x | dense-matrix | no | 4.01e-15 |
| 33 | 8320 | complex | cosine | inverse | 0.3695 | 1.5195 | 0.5730 | 0.645x | dense-matrix | no | 5.33e-15 |
| 33 | 8320 | complex | sine | forward | 0.3051 | 1.3831 | 0.2562 | 1.191x | fftw-allocating | yes | 3.98e-15 |
| 33 | 8320 | complex | sine | inverse | 0.3567 | 1.3860 | 0.6028 | 0.592x | dense-matrix | no | 5.25e-15 |
| 33 | 33024 | real | cosine | forward | 0.2963 | 3.6668 | 1.2260 | 0.242x | dense-matrix | no | 4.91e-15 |
| 33 | 33024 | real | cosine | inverse | 0.2794 | 4.3986 | 1.4392 | 0.194x | dense-matrix | no | 4.78e-15 |
| 33 | 33024 | real | sine | forward | 0.2457 | 4.2095 | 0.5075 | 0.484x | dense-matrix | no | 4.87e-15 |
| 33 | 33024 | real | sine | inverse | 0.2843 | 4.1457 | 1.4740 | 0.193x | dense-matrix | no | 4.63e-15 |
| 33 | 33024 | complex | cosine | forward | 3.1563 | 5.0135 | 1.9673 | 1.604x | fftw-allocating | yes | 3.64e-15 |
| 33 | 33024 | complex | cosine | inverse | 3.2012 | 5.1450 | 2.2317 | 1.434x | fftw-allocating | yes | 4.79e-15 |
| 33 | 33024 | complex | sine | forward | 2.7343 | 4.3487 | 1.0869 | 2.516x | fftw-allocating | yes | 5.02e-15 |
| 33 | 33024 | complex | sine | inverse | 1.3140 | 4.3393 | 2.0087 | 0.654x | dense-matrix | no | 5.39e-15 |
| 33 | 131584 | real | cosine | forward | 3.7408 | 13.4043 | 4.7609 | 0.786x | dense-matrix | no | 4.83e-15 |
| 33 | 131584 | real | cosine | inverse | 3.8030 | 16.2172 | 5.9715 | 0.637x | dense-matrix | no | 4.25e-15 |
| 33 | 131584 | real | sine | forward | 3.2955 | 15.8746 | 1.4364 | 2.294x | fftw-allocating | yes | 4.83e-15 |
| 33 | 131584 | real | sine | inverse | 3.6853 | 15.2611 | 5.7331 | 0.643x | dense-matrix | no | 5.39e-15 |
| 33 | 131584 | complex | cosine | forward | 12.8015 | 18.4962 | 6.6756 | 1.918x | fftw-allocating | yes | 4.21e-15 |
| 33 | 131584 | complex | cosine | inverse | 12.6702 | 20.9298 | 10.8809 | 1.164x | fftw-allocating | yes | 4.75e-15 |
| 33 | 131584 | complex | sine | forward | 11.1294 | 16.4020 | 3.3722 | 3.300x | fftw-allocating | yes | 4.59e-15 |
| 33 | 131584 | complex | sine | inverse | 12.4365 | 19.0879 | 7.7305 | 1.609x | fftw-allocating | yes | 5.15e-15 |
| 65 | 1 | real | cosine | forward | 0.0019 | 0.0898 | 0.1620 | 0.012x | dense-matrix | no | 1.01e-14 |
| 65 | 1 | real | cosine | inverse | 0.0012 | 0.1190 | 0.1665 | 0.007x | dense-matrix | no | 5.66e-15 |
| 65 | 1 | real | sine | forward | 0.0011 | 0.0944 | 0.1455 | 0.008x | dense-matrix | no | 5.43e-15 |
| 65 | 1 | real | sine | inverse | 0.0010 | 0.0720 | 0.1370 | 0.007x | dense-matrix | no | 5.81e-15 |
| 65 | 1 | complex | cosine | forward | 0.0129 | 0.0631 | 0.0926 | 0.139x | dense-matrix | no | 4.89e-15 |
| 65 | 1 | complex | cosine | inverse | 0.0047 | 0.0735 | 0.0830 | 0.057x | dense-matrix | no | 6.36e-15 |
| 65 | 1 | complex | sine | forward | 0.0038 | 0.0529 | 0.1441 | 0.027x | dense-matrix | no | 6.74e-15 |
| 65 | 1 | complex | sine | inverse | 0.0035 | 0.0488 | 0.1542 | 0.023x | dense-matrix | no | 7.1e-15 |
| 65 | 8320 | real | cosine | forward | 0.2494 | 1.8708 | 0.7175 | 0.348x | dense-matrix | no | 7.57e-15 |
| 65 | 8320 | real | cosine | inverse | 0.2535 | 2.3930 | 0.8007 | 0.317x | dense-matrix | no | 9.66e-15 |
| 65 | 8320 | real | sine | forward | 0.2151 | 2.1600 | 0.2675 | 0.804x | dense-matrix | no | 8.01e-15 |
| 65 | 8320 | real | sine | inverse | 0.2341 | 2.0277 | 0.7651 | 0.306x | dense-matrix | no | 8.72e-15 |
| 65 | 8320 | complex | cosine | forward | 0.9058 | 2.5505 | 1.0340 | 0.876x | dense-matrix | no | 7.81e-15 |
| 65 | 8320 | complex | cosine | inverse | 0.9106 | 3.3872 | 1.1057 | 0.824x | dense-matrix | no | 1e-14 |
| 65 | 8320 | complex | sine | forward | 0.8225 | 2.4475 | 0.5989 | 1.373x | fftw-allocating | yes | 9.08e-15 |
| 65 | 8320 | complex | sine | inverse | 0.8767 | 2.7854 | 0.9937 | 0.882x | dense-matrix | no | 8.19e-15 |
| 65 | 33024 | real | cosine | forward | 1.7920 | 6.5088 | 2.6980 | 0.664x | dense-matrix | no | 7.85e-15 |
| 65 | 33024 | real | cosine | inverse | 1.7985 | 7.6356 | 3.0568 | 0.588x | dense-matrix | no | 1.09e-14 |
| 65 | 33024 | real | sine | forward | 1.6790 | 7.1175 | 1.2040 | 1.394x | fftw-allocating | yes | 7.75e-15 |
| 65 | 33024 | real | sine | inverse | 0.8800 | 7.1819 | 2.8528 | 0.308x | dense-matrix | no | 9.15e-15 |
| 65 | 33024 | complex | cosine | forward | 6.3228 | 9.7682 | 3.5905 | 1.761x | fftw-allocating | yes | 9.04e-15 |
| 65 | 33024 | complex | cosine | inverse | 6.3876 | 11.9657 | 4.5823 | 1.394x | fftw-allocating | yes | 1.04e-14 |
| 65 | 33024 | complex | sine | forward | 5.8781 | 8.8456 | 2.1748 | 2.703x | fftw-allocating | yes | 7.72e-15 |
| 65 | 33024 | complex | sine | inverse | 6.3131 | 10.8291 | 3.7954 | 1.663x | fftw-allocating | yes | 9.53e-15 |
| 65 | 131584 | real | cosine | forward | 7.1040 | 25.6851 | 10.1940 | 0.697x | dense-matrix | no | 8.23e-15 |
| 65 | 131584 | real | cosine | inverse | 7.0556 | 29.8659 | 12.8150 | 0.551x | dense-matrix | no | 9.56e-15 |
| 65 | 131584 | real | sine | forward | 6.6015 | 26.7592 | 3.4901 | 1.891x | fftw-allocating | yes | 8.97e-15 |
| 65 | 131584 | real | sine | inverse | 7.0048 | 26.8702 | 11.3098 | 0.619x | dense-matrix | no | 9.72e-15 |
| 65 | 131584 | complex | cosine | forward | 25.3928 | 36.2951 | 14.2963 | 1.776x | fftw-allocating | yes | 7.39e-15 |
| 65 | 131584 | complex | cosine | inverse | 25.4872 | 41.7585 | 20.7547 | 1.228x | fftw-allocating | yes | 9.13e-15 |
| 65 | 131584 | complex | sine | forward | 23.4845 | 34.7080 | 7.6744 | 3.060x | fftw-allocating | yes | 8.67e-15 |
| 65 | 131584 | complex | sine | inverse | 25.1947 | 35.1296 | 15.2597 | 1.651x | fftw-allocating | yes | 8.27e-15 |
| 129 | 1 | real | cosine | forward | 0.0013 | 0.0292 | 0.1394 | 0.009x | dense-matrix | no | 1.35e-14 |
| 129 | 1 | real | cosine | inverse | 0.0015 | 0.0783 | 0.1429 | 0.010x | dense-matrix | no | 1.11e-14 |
| 129 | 1 | real | sine | forward | 0.0014 | 0.0227 | 0.1300 | 0.011x | dense-matrix | no | 1.14e-14 |
| 129 | 1 | real | sine | inverse | 0.0015 | 0.0683 | 0.1349 | 0.011x | dense-matrix | no | 1.22e-14 |
| 129 | 1 | complex | cosine | forward | 0.0079 | 0.0537 | 0.0673 | 0.118x | dense-matrix | no | 1.6e-14 |
| 129 | 1 | complex | cosine | inverse | 0.0077 | 0.0430 | 0.0664 | 0.116x | dense-matrix | no | 1.23e-14 |
| 129 | 1 | complex | sine | forward | 0.0075 | 0.0428 | 0.1510 | 0.050x | dense-matrix | no | 9.02e-15 |
| 129 | 1 | complex | sine | inverse | 0.0072 | 0.0548 | 0.1489 | 0.048x | dense-matrix | no | 1.45e-14 |
| 129 | 8320 | real | cosine | forward | 0.7201 | 3.8564 | 1.3221 | 0.545x | dense-matrix | no | 1.25e-14 |
| 129 | 8320 | real | cosine | inverse | 0.7223 | 4.5687 | 1.4679 | 0.492x | dense-matrix | no | 1.98e-14 |
| 129 | 8320 | real | sine | forward | 0.6630 | 4.2736 | 0.5008 | 1.324x | fftw-allocating | yes | 1.77e-14 |
| 129 | 8320 | real | sine | inverse | 0.7039 | 4.1888 | 1.3783 | 0.511x | dense-matrix | no | 1.77e-14 |
| 129 | 8320 | complex | cosine | forward | 4.0053 | 4.9270 | 1.8730 | 2.138x | fftw-allocating | yes | 1.41e-14 |
| 129 | 8320 | complex | cosine | inverse | 3.9953 | 6.5826 | 2.0045 | 1.993x | fftw-allocating | yes | 1.77e-14 |
| 129 | 8320 | complex | sine | forward | 3.8596 | 5.0383 | 1.1375 | 3.393x | fftw-allocating | yes | 1.69e-14 |
| 129 | 8320 | complex | sine | inverse | 3.9396 | 5.4588 | 1.7536 | 2.247x | fftw-allocating | yes | 1.68e-14 |
| 129 | 33024 | real | cosine | forward | 4.0862 | 13.3993 | 4.8907 | 0.836x | dense-matrix | no | 1.58e-14 |
| 129 | 33024 | real | cosine | inverse | 4.1027 | 15.8313 | 5.8300 | 0.704x | dense-matrix | no | 1.87e-14 |
| 129 | 33024 | real | sine | forward | 3.9120 | 15.1738 | 1.6410 | 2.384x | fftw-allocating | yes | 1.78e-14 |
| 129 | 33024 | real | sine | inverse | 4.0550 | 13.9747 | 5.0882 | 0.797x | dense-matrix | no | 1.53e-14 |
| 129 | 33024 | complex | cosine | forward | 15.8242 | 20.5903 | 6.3227 | 2.503x | fftw-allocating | yes | 1.75e-14 |
| 129 | 33024 | complex | cosine | inverse | 15.8271 | 25.2731 | 7.4743 | 2.118x | fftw-allocating | yes | 1.76e-14 |
| 129 | 33024 | complex | sine | forward | 15.2343 | 20.4847 | 3.0360 | 5.018x | fftw-allocating | yes | 1.69e-14 |
| 129 | 33024 | complex | sine | inverse | 15.6843 | 22.6089 | 6.3859 | 2.456x | fftw-allocating | yes | 1.9e-14 |
| 129 | 131584 | real | cosine | forward | 16.4627 | 51.3842 | 19.0014 | 0.866x | dense-matrix | no | 1.59e-14 |
| 129 | 131584 | real | cosine | inverse | 16.3768 | 59.4942 | 23.8860 | 0.686x | dense-matrix | no | 1.72e-14 |
| 129 | 131584 | real | sine | forward | 15.5287 | 56.9429 | 5.7557 | 2.698x | fftw-allocating | yes | 1.72e-14 |
| 129 | 131584 | real | sine | inverse | 16.1440 | 51.3697 | 20.3423 | 0.794x | dense-matrix | no | 1.7e-14 |
| 129 | 131584 | complex | cosine | forward | 62.8190 | 76.3305 | 23.2488 | 2.702x | fftw-allocating | yes | 1.52e-14 |
| 129 | 131584 | complex | cosine | inverse | 62.6608 | 95.2333 | 30.8585 | 2.031x | fftw-allocating | yes | 1.8e-14 |
| 129 | 131584 | complex | sine | forward | 60.5328 | 78.3901 | 11.6130 | 5.213x | fftw-allocating | yes | 1.78e-14 |
| 129 | 131584 | complex | sine | inverse | 61.9920 | 82.0303 | 26.0927 | 2.376x | fftw-allocating | yes | 1.82e-14 |
| 257 | 1 | real | cosine | forward | 0.0025 | 0.0212 | 0.1479 | 0.017x | dense-matrix | no | 2.86e-14 |
| 257 | 1 | real | cosine | inverse | 0.0025 | 0.0182 | 0.1497 | 0.017x | dense-matrix | no | 1.88e-14 |
| 257 | 1 | real | sine | forward | 0.0020 | 0.0242 | 0.1372 | 0.015x | dense-matrix | no | 3.19e-14 |
| 257 | 1 | real | sine | inverse | 0.0020 | 0.0693 | 0.1359 | 0.015x | dense-matrix | no | 2.85e-14 |
| 257 | 1 | complex | cosine | forward | 0.0255 | 0.0364 | 0.0640 | 0.399x | dense-matrix | no | 3.11e-14 |
| 257 | 1 | complex | cosine | inverse | 0.0229 | 0.0336 | 0.0587 | 0.390x | dense-matrix | no | 2.73e-14 |
| 257 | 1 | complex | sine | forward | 0.0238 | 0.0198 | 0.1434 | 0.138x | fft-extension | no | 1.88e-14 |
| 257 | 1 | complex | sine | inverse | 0.0227 | 0.0191 | 0.1483 | 0.129x | fft-extension | no | 2.84e-14 |
| 257 | 8320 | real | cosine | forward | 2.9972 | 8.0236 | 2.4986 | 1.200x | fftw-allocating | yes | 3.41e-14 |
| 257 | 8320 | real | cosine | inverse | 2.9907 | 9.2975 | 2.8702 | 1.042x | fftw-allocating | no | 3.38e-14 |
| 257 | 8320 | real | sine | forward | 2.9087 | 8.6754 | 1.0426 | 2.790x | fftw-allocating | yes | 3.28e-14 |
| 257 | 8320 | real | sine | inverse | 2.9971 | 8.4070 | 2.5516 | 1.175x | fftw-allocating | yes | 3.22e-14 |
| 257 | 8320 | complex | cosine | forward | 12.1715 | 9.8477 | 3.3505 | 2.939x | fftw-allocating | yes | 3.06e-14 |
| 257 | 8320 | complex | cosine | inverse | 12.4237 | 13.7399 | 3.8465 | 3.230x | fftw-allocating | yes | 3.62e-14 |
| 257 | 8320 | complex | sine | forward | 11.9905 | 10.7070 | 1.7460 | 6.132x | fftw-allocating | yes | 3.55e-14 |
| 257 | 8320 | complex | sine | inverse | 12.1607 | 11.2575 | 3.1876 | 3.532x | fftw-allocating | yes | 3.26e-14 |
| 257 | 33024 | real | cosine | forward | 11.8962 | 31.6327 | 9.7315 | 1.222x | fftw-allocating | yes | 3.29e-14 |
| 257 | 33024 | real | cosine | inverse | 11.8949 | 35.7627 | 11.2023 | 1.062x | fftw-allocating | no | 3.5e-14 |
| 257 | 33024 | real | sine | forward | 11.5761 | 33.7662 | 3.0235 | 3.829x | fftw-allocating | yes | 3.59e-14 |
| 257 | 33024 | real | sine | inverse | 11.8491 | 30.9119 | 9.7614 | 1.214x | fftw-allocating | yes | 3.4e-14 |
| 257 | 33024 | complex | cosine | forward | 47.3010 | 40.2090 | 12.1877 | 3.299x | fftw-allocating | yes | 3.06e-14 |
| 257 | 33024 | complex | cosine | inverse | 47.3030 | 50.6475 | 14.5324 | 3.255x | fftw-allocating | yes | 3.49e-14 |
| 257 | 33024 | complex | sine | forward | 46.5652 | 40.8020 | 5.9520 | 6.855x | fftw-allocating | yes | 3.71e-14 |
| 257 | 33024 | complex | sine | inverse | 47.3324 | 42.5253 | 12.3773 | 3.436x | fftw-allocating | yes | 3.49e-14 |
| 257 | 131584 | real | cosine | forward | 47.2431 | 121.2950 | 38.5484 | 1.226x | fftw-allocating | yes | 2.37e-14 |
| 257 | 131584 | real | cosine | inverse | 47.1585 | 140.9422 | 45.0524 | 1.047x | fftw-allocating | no | 3.37e-14 |
| 257 | 131584 | real | sine | forward | 45.8560 | 131.8640 | 11.3925 | 4.025x | fftw-allocating | yes | 3.36e-14 |
| 257 | 131584 | real | sine | inverse | 46.8755 | 119.1037 | 39.6837 | 1.181x | fftw-allocating | yes | 3.53e-14 |
| 257 | 131584 | complex | cosine | forward | 187.5161 | 152.3312 | 46.9271 | 3.246x | fftw-allocating | yes | 2.84e-14 |
| 257 | 131584 | complex | cosine | inverse | 187.8315 | 194.2309 | 56.8021 | 3.307x | fftw-allocating | yes | 3.49e-14 |
| 257 | 131584 | complex | sine | forward | 184.5517 | 157.0677 | 22.3410 | 7.030x | fftw-allocating | yes | 3.34e-14 |
| 257 | 131584 | complex | sine | inverse | 187.2890 | 165.6787 | 50.3469 | 3.291x | fftw-allocating | yes | 3.56e-14 |
| 513 | 1 | real | cosine | forward | 0.0053 | 0.0188 | 0.4004 | 0.013x | dense-matrix | no | 5.18e-14 |
| 513 | 1 | real | cosine | inverse | 0.0058 | 0.0188 | 0.3850 | 0.015x | dense-matrix | no | 6.21e-14 |
| 513 | 1 | real | sine | forward | 0.0056 | 0.0250 | 0.3984 | 0.014x | dense-matrix | no | 6.6e-14 |
| 513 | 1 | real | sine | inverse | 0.0063 | 0.0258 | 0.3840 | 0.016x | dense-matrix | no | 5.49e-14 |
| 513 | 1 | complex | cosine | forward | 0.0807 | 0.0222 | 0.3302 | 0.067x | fft-extension | no | 5.09e-14 |
| 513 | 1 | complex | cosine | inverse | 0.0825 | 0.0302 | 0.3352 | 0.090x | fft-extension | no | 4.29e-14 |
| 513 | 1 | complex | sine | forward | 0.0821 | 0.0226 | 0.4178 | 0.054x | fft-extension | no | 5.95e-14 |
| 513 | 1 | complex | sine | inverse | 0.0811 | 0.0210 | 0.4103 | 0.051x | fft-extension | no | 5.96e-14 |
| 513 | 8320 | real | cosine | forward | 10.1387 | 16.0352 | 5.1495 | 1.969x | fftw-allocating | yes | 6.96e-14 |
| 513 | 8320 | real | cosine | inverse | 10.0604 | 18.1669 | 5.6050 | 1.795x | fftw-allocating | yes | 6.63e-14 |
| 513 | 8320 | real | sine | forward | 10.1463 | 17.2132 | 1.8219 | 5.569x | fftw-allocating | yes | 6.94e-14 |
| 513 | 8320 | real | sine | inverse | 10.2814 | 15.8341 | 5.0361 | 2.042x | fftw-allocating | yes | 8.26e-14 |
| 513 | 8320 | complex | cosine | forward | 44.5924 | 17.7213 | 6.2650 | 2.829x | fftw-allocating | yes | 6.71e-14 |
| 513 | 8320 | complex | cosine | inverse | 44.9318 | 23.8872 | 7.0470 | 3.390x | fftw-allocating | yes | 7.31e-14 |
| 513 | 8320 | complex | sine | forward | 44.2013 | 18.3050 | 3.4542 | 5.299x | fftw-allocating | yes | 6.94e-14 |
| 513 | 8320 | complex | sine | inverse | 39.5418 | 19.1776 | 6.5106 | 2.946x | fftw-allocating | yes | 6.69e-14 |
| 513 | 33024 | real | cosine | forward | 39.5481 | 61.0900 | 19.4890 | 2.029x | fftw-allocating | yes | 6.55e-14 |
| 513 | 33024 | real | cosine | inverse | 39.4933 | 70.9726 | 21.8346 | 1.809x | fftw-allocating | yes | 6.33e-14 |
| 513 | 33024 | real | sine | forward | 38.9738 | 67.7533 | 5.7363 | 6.794x | fftw-allocating | yes | 6.63e-14 |
| 513 | 33024 | real | sine | inverse | 39.3564 | 61.4258 | 19.5232 | 2.016x | fftw-allocating | yes | 6.95e-14 |
| 513 | 33024 | complex | cosine | forward | 176.3744 | 66.4975 | 23.2937 | 2.855x | fftw-allocating | yes | 5.5e-14 |
| 513 | 33024 | complex | cosine | inverse | 176.4038 | 85.0620 | 26.7707 | 3.177x | fftw-allocating | yes | 6.49e-14 |
| 513 | 33024 | complex | sine | forward | 174.7693 | 67.2342 | 11.2178 | 5.994x | fftw-allocating | yes | 6.65e-14 |
| 513 | 33024 | complex | sine | inverse | 155.8582 | 74.9755 | 24.4093 | 3.072x | fftw-allocating | yes | 6.19e-14 |
| 513 | 131584 | real | cosine | forward | 157.1273 | 234.4391 | 77.2519 | 2.034x | fftw-allocating | yes | 6.11e-14 |
| 513 | 131584 | real | cosine | inverse | 157.7027 | 268.7749 | 86.3729 | 1.826x | fftw-allocating | yes | 6.71e-14 |
| 513 | 131584 | real | sine | forward | 154.8502 | 263.4314 | 21.7963 | 7.104x | fftw-allocating | yes | 6.64e-14 |
| 513 | 131584 | real | sine | inverse | 162.6783 | 247.0361 | 82.2742 | 1.977x | fftw-allocating | yes | 6.48e-14 |
| 513 | 131584 | complex | cosine | forward | 737.8106 | 286.8159 | 94.0057 | 3.051x | fftw-allocating | yes | 6.24e-14 |
| 513 | 131584 | complex | cosine | inverse | 711.4250 | 335.6863 | 104.8187 | 3.203x | fftw-allocating | yes | 6.73e-14 |
| 513 | 131584 | complex | sine | forward | 700.7991 | 286.2197 | 41.5570 | 6.887x | fftw-allocating | yes | 7.41e-14 |
| 513 | 131584 | complex | sine | inverse | 621.7658 | 290.0624 | 102.8427 | 2.820x | fftw-allocating | yes | 7.37e-14 |

## Bounded eligibility

| Nz | Type | Transform | Direction | Eligible batch intervals |
|---:|---|---|---|---|
| 33 | real | cosine | forward | none |
| 33 | real | cosine | inverse | none |
| 33 | real | sine | forward | 131584-131584 |
| 33 | real | sine | inverse | none |
| 33 | complex | cosine | forward | 33024-131584 |
| 33 | complex | cosine | inverse | 33024-131584 |
| 33 | complex | sine | forward | 8320-131584 |
| 33 | complex | sine | inverse | 131584-131584 |
| 65 | real | cosine | forward | none |
| 65 | real | cosine | inverse | none |
| 65 | real | sine | forward | 33024-131584 |
| 65 | real | sine | inverse | none |
| 65 | complex | cosine | forward | 33024-131584 |
| 65 | complex | cosine | inverse | 33024-131584 |
| 65 | complex | sine | forward | 8320-131584 |
| 65 | complex | sine | inverse | 33024-131584 |
| 129 | real | cosine | forward | none |
| 129 | real | cosine | inverse | none |
| 129 | real | sine | forward | 8320-131584 |
| 129 | real | sine | inverse | none |
| 129 | complex | cosine | forward | 8320-131584 |
| 129 | complex | cosine | inverse | 8320-131584 |
| 129 | complex | sine | forward | 8320-131584 |
| 129 | complex | sine | inverse | 8320-131584 |
| 257 | real | cosine | forward | 8320-131584 |
| 257 | real | cosine | inverse | none |
| 257 | real | sine | forward | 8320-131584 |
| 257 | real | sine | inverse | 8320-131584 |
| 257 | complex | cosine | forward | 8320-131584 |
| 257 | complex | cosine | inverse | 8320-131584 |
| 257 | complex | sine | forward | 8320-131584 |
| 257 | complex | sine | inverse | 8320-131584 |
| 513 | real | cosine | forward | 8320-131584 |
| 513 | real | cosine | inverse | 8320-131584 |
| 513 | real | sine | forward | 8320-131584 |
| 513 | real | sine | inverse | 8320-131584 |
| 513 | complex | cosine | forward | 8320-131584 |
| 513 | complex | cosine | inverse | 8320-131584 |
| 513 | complex | sine | forward | 8320-131584 |
| 513 | complex | sine | inverse | 8320-131584 |

## Timing boundaries

- Complete-call time includes the MATLAB method or reference transform call.
- FFTW kernel time contains fftw_execute_r2r only.
- FFTW pipeline time contains kernel execution plus GL normalization and endpoint handling.
- Input generation, matrix construction, planning, preallocation, and diagnostic retrieval are excluded.
