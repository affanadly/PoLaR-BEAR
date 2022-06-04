# Python Language Radio Burst Emission Automatic Roger (PoLaR BEAR)
An FRB detection code written for my undergraduate research project. This Python code is written based on the Burst Emission Automatic Roger (BEAR) written in `C++` (Men et al., 2019). An explanation of the data analysis process is available in the same article as well as in my [thesis](https://github.com/affanadly/PoLaR-BEAR/blob/main/thesis/Research_Thesis__Fast_Radio_Bursts.pdf). 

# Notes
* This repository is only intended as an archive, no further development will be done. This code was finalized in June 2021.
* PoLaR BEAR was written for me to learn and understand FRB detection techniques. Therefore, it is mostly inefficient.
* PoLaR BEAR was written as a functional code rather than object oriented. 
* Inputs and analysis parameters are not taken using arguments, the code has to be modified to adjust them. 
* New filterbank header items (which are added after June 2021) which are not available in the header dictionary will cause the code to fail. Add the new attributes to get it to work again. 
* Fake FRB filterbank data can be generated with using [FRBFakeRandom.py](https://github.com/affanadly/PoLaR-BEAR/blob/main/FRBFakeRandom.py).

# Example output
The data used here is the FRB010124 (source: [Parkes Archival FRB Data](https://data-portal.hpc.swin.edu.au/dataset/parkes-frbs-archival-data)).

* Main plot
![Main plot](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_Plot.jpeg)
* Candidate pulse profiles
![Candidate pulse profiles](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_Candidates.jpeg)
* Individual candidate plots
![Individual candidate plot 1](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_C1.jpeg)
![Individual candidate plot 2](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_C2.jpeg)
![Individual candidate plot 3](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_C3.jpeg)
![Individual candidate plot 4](https://github.com/affanadly/PoLaR-BEAR/blob/main/example/FRB_51934.019976851850515_C4.jpeg)
* Parameter output
```
Source: Unknown
Time (MJD): 51934.019976851850515

Downsampling coefficient: 20
Start: 0.00000000
End: 0.20000000
Threshold: 41.99597534

Zapped frequencies:
1435.00 - 1440.00 MHz
1495.00 - 1505.00 MHz

 Candidate    DM (cm-3 pc)    Width (s)    Time (s)      S       SNR
-----------  --------------  -----------  ----------  -------  -------
     1            790          0.0075      28.8025    507.762  22.5336
     2             1           0.0025       22.535    95.7433  9.78485
     3            671          0.1225       28.855    71.634   8.46369
     4            981          0.1225       28.735    46.0027  6.78253
```

# Reference
* Men, Y. P., Luo, R., Chen, M. Z., Hao, L. F., Lee, K. J., Li, J., Li, Z. X., Liu, Z. Y., Pei, X., Wen, Z. G., Wu, J. J., Xu, Y. H., Xu, R. X., Yuan, J. P., & Zhang, C. F. (2019). Piggyback search for fast radio bursts using Nanshan 26 m and Kunming 40 m radio telescopes – I. Observing and data analysis systems, discovery of a mysterious peryton. In Monthly Notices of the Royal Astronomical Society (Vol. 488, Issue 3, pp. 3957–3971). Oxford University Press (OUP). https://doi.org/10.1093/mnras/stz1931
