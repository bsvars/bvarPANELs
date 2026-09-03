# Provides posterior estimation summary for Bayesian Hierarchical Panel Vector Autoregressions with group-specific global prior

Provides posterior mean, standard deviations, as well as 5 and 95
percentiles of the parameters for all `C` countries.

## Usage

``` r
# S3 method for class 'PosteriorBVARGROUPPRIORPANEL'
summary(object, ...)
```

## Arguments

- object:

  an object of class `PosteriorBVARGROUPPRIORPANEL` obtained using the
  [`estimate()`](https://bsvars.org/bsvars/reference/estimate.html)
  function applied to Vector Autoregressions containing draws from the
  posterior distribution of the parameters.

- ...:

  additional arguments affecting the summary produced.

## Value

A list reporting the posterior mean, standard deviations, as well as 5
and 95 percentiles of the country-specific parameters.

## See also

[`estimate.BVARGROUPPRIORPANEL`](http://bsvars.org/bpvars/dev/reference/estimate.BVARGROUPPRIORPANEL.md),
[`specify_bvarGroupPriorPANEL`](http://bsvars.org/bpvars/dev/reference/specify_bvarGroupPriorPANEL.md)

## Author

Tomasz Woźniak <wozniak.tom@pm.me>

## Examples

``` r
# specify the model
specification = specify_bvarGroupPriorPANEL$new(
      data = ilo_dynamic_panel[1:5],
      group_allocation = country_grouping_region[1:5]
)
#> Country groupings have been pre-specified and will not be estimated.
burn_in       = estimate(specification, 5)             # run the burn-in
#> **************************************************|
#>  bpvars: Forecasting with Bayesian Panel VARs     |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
posterior     = estimate(burn_in, 5)                   # estimate the model
#> **************************************************|
#>  bpvars: Forecasting with Bayesian Panel VARs     |
#> **************************************************|
#>  Progress of the MCMC simulation for 5 draws
#>     Every draw is saved via MCMC thinning
#>  Press Esc to interrupt the computations
#> **************************************************|
summary(posterior)
#> $AFG
#> $AFG$Sigma
#> $AFG$Sigma$equation1
#>                   mean          sd 5% quantile 95% quantile
#> Sigma[1,1] 0.004657898 0.001248489 0.003394889  0.006163106
#> 
#> $AFG$Sigma$equation2
#>                    mean          sd  5% quantile  95% quantile
#> Sigma[2,1] -0.004989466 0.004395443 -0.009705882 -0.0003755101
#> Sigma[2,2]  0.186980638 0.031041126  0.149212377  0.2184459294
#> 
#> $AFG$Sigma$equation3
#>                     mean          sd 5% quantile 95% quantile
#> Sigma[3,1] -0.0004865243 0.006782181 -0.00548591  0.008713639
#> Sigma[3,2] -0.1625804381 0.060772303 -0.22527925 -0.093873436
#> Sigma[3,3]  0.3464548142 0.074112040  0.27515830  0.439146735
#> 
#> $AFG$Sigma$equation4
#>                    mean          sd  5% quantile 95% quantile
#> Sigma[4,1] -0.003868859 0.006111286 -0.009253865  0.004277054
#> Sigma[4,2] -0.092144419 0.055379137 -0.146391328 -0.032189025
#> Sigma[4,3]  0.310139208 0.054491037  0.260735549  0.381102943
#> Sigma[4,4]  0.306034174 0.039350227  0.280594976  0.359521603
#> 
#> 
#> $AFG$A
#> $AFG$A$equation1
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1  0.9771281 0.02569019   0.9505739    1.0083647
#> lag1_var2  0.2107912 0.01906929   0.1862454    0.2274251
#> lag1_var3  0.4991445 0.01798942   0.4837799    0.5224646
#> lag1_var4 -0.4374319 0.01791093  -0.4589854   -0.4192865
#> const     -2.1076326 0.70574431  -2.9529858   -1.4540891
#> 
#> $AFG$A$equation2
#>                 mean         sd  5% quantile 95% quantile
#> lag1_var1  0.1925272 0.15881830  0.004157616    0.3695079
#> lag1_var2  0.7792215 0.09200682  0.653188505    0.8481498
#> lag1_var3 -0.5321093 0.13636239 -0.705915183   -0.4071184
#> lag1_var4  0.4477190 0.12232281  0.337353539    0.5945298
#> const     -0.5679611 4.53951618 -4.878802097    5.3323664
#> 
#> $AFG$A$equation3
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1  0.0972481 0.07282411  0.04598929   0.19725349
#> lag1_var2 -0.2962999 0.11079621 -0.41528935  -0.15991764
#> lag1_var3  1.1940223 0.13324719  1.01409163   1.30198856
#> lag1_var4 -0.2706983 0.13444994 -0.38726531  -0.09062673
#> const      4.2649269 3.57929430 -0.64698843   6.87408798
#> 
#> $AFG$A$equation4
#>                  mean         sd 5% quantile 95% quantile
#> lag1_var1  0.19837625 0.07713060  0.10114100   0.27089191
#> lag1_var2  0.08695486 0.08887674 -0.02632473   0.17585699
#> lag1_var3  1.09211281 0.10821829  0.95732638   1.19725844
#> lag1_var4 -0.06296876 0.12057217 -0.19647334   0.07782666
#> const     -2.66305782 2.60443657 -5.64516029   0.15591057
#> 
#> 
#> 
#> $AGO
#> $AGO$Sigma
#> $AGO$Sigma$equation1
#>                  mean          sd 5% quantile 95% quantile
#> Sigma[1,1] 0.01421447 0.005805898 0.008228763   0.02137606
#> 
#> $AGO$Sigma$equation2
#>                    mean         sd 5% quantile 95% quantile
#> Sigma[2,1] -0.007348488 0.02595655 -0.03881812   0.02117339
#> Sigma[2,2]  0.255900225 0.11344521  0.14620344   0.40254368
#> 
#> $AGO$Sigma$equation3
#>                    mean         sd 5% quantile 95% quantile
#> Sigma[3,1] -0.003858301 0.02149805  -0.0275049    0.0208399
#> Sigma[3,2] -0.219508805 0.12126356  -0.3775213   -0.1267110
#> Sigma[3,3]  0.450783391 0.18692119   0.2759204    0.6730613
#> 
#> $AGO$Sigma$equation4
#>                   mean          sd 5% quantile 95% quantile
#> Sigma[4,1] -0.01325469 0.005530217 -0.01758409 -0.005848214
#> Sigma[4,2] -0.02714567 0.093056173 -0.12889329  0.084810773
#> Sigma[4,3]  0.34322337 0.144461646  0.19489846  0.526851977
#> Sigma[4,4]  0.39663153 0.103748472  0.31095876  0.535609299
#> 
#> 
#> $AGO$A
#> $AGO$A$equation1
#>                 mean        sd 5% quantile 95% quantile
#> lag1_var1  1.0728659 0.0291689   1.0477943    1.1114669
#> lag1_var2 -0.5219895 0.1649731  -0.6886576   -0.3118732
#> lag1_var3 -0.5948180 0.2221713  -0.8281364   -0.3186681
#> lag1_var4  0.5911684 0.2031935   0.3360705    0.8033620
#> const     -0.4257331 1.6825932  -1.8653155    1.8155613
#> 
#> $AGO$A$equation2
#>                 mean        sd 5% quantile 95% quantile
#> lag1_var1 -0.2187796 0.2461240  -0.4402469    0.0938496
#> lag1_var2  1.0748833 0.8186163   0.3701250    2.1721424
#> lag1_var3  0.5322424 1.0419790  -0.2674114    1.9495803
#> lag1_var4 -0.3466586 0.9970829  -1.6842554    0.4241860
#> const     -3.4408983 3.0749807  -7.4582747   -0.5473634
#> 
#> $AGO$A$equation3
#>                  mean        sd 5% quantile 95% quantile
#> lag1_var1 -0.08796298 0.2857986  -0.4103201    0.2553907
#> lag1_var2  1.81686140 0.9460023   0.6996367    2.8304726
#> lag1_var3  2.76852733 1.1322981   1.3359427    3.8446927
#> lag1_var4 -1.93661988 1.0815348  -2.9817225   -0.5573837
#> const      7.50660322 4.3679378   1.7821225   11.5962092
#> 
#> $AGO$A$equation4
#>                 mean        sd 5% quantile 95% quantile
#> lag1_var1 -0.3613178 0.1796619  -0.4795812   -0.1149353
#> lag1_var2  2.8676831 0.5229266   2.2218996    3.3836675
#> lag1_var3  3.4422590 0.6256053   2.5898220    3.9115371
#> lag1_var4 -2.3779603 0.5434040  -2.8507909   -1.6678264
#> const      0.1828613 4.7944138  -5.3804846    5.7751857
#> 
#> 
#> 
#> $ALB
#> $ALB$Sigma
#> $ALB$Sigma$equation1
#>                   mean          sd 5% quantile 95% quantile
#> Sigma[1,1] 0.005501371 0.002343237  0.00248115  0.007627412
#> 
#> $ALB$Sigma$equation2
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[2,1] -0.03953674 0.02450947 -0.06480436  -0.01328715
#> Sigma[2,2]  6.74397657 1.42303653  5.18371359   8.32665223
#> 
#> $ALB$Sigma$equation3
#>                    mean         sd 5% quantile 95% quantile
#> Sigma[3,1] -0.006885861 0.03605658 -0.05570963   0.01370019
#> Sigma[3,2] -3.815953767 1.36860856 -5.45353852  -2.43679581
#> Sigma[3,3]  4.997698834 0.67200979  4.33891485   5.83440188
#> 
#> $ALB$Sigma$equation4
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[4,1] -0.03774451 0.04721912  -0.0990949  0.003057428
#> Sigma[4,2]  0.81882032 0.79565294  -0.1787429  1.660246369
#> Sigma[4,3]  2.90488217 1.37501053   1.8019527  4.745199207
#> Sigma[4,4]  4.12816913 2.00052688   2.3634579  6.561301223
#> 
#> 
#> $ALB$A
#> $ALB$A$equation1
#>                   mean         sd 5% quantile 95% quantile
#> lag1_var1  1.054818152 0.03023398  1.02772709   1.09136835
#> lag1_var2  0.014522601 0.04336567 -0.04397189   0.04888939
#> lag1_var3  0.008931424 0.07238784 -0.08871198   0.06677327
#> lag1_var4 -0.008586858 0.05992156 -0.05548270   0.07230708
#> const     -1.375338817 0.25828911 -1.65848681  -1.08618953
#> 
#> $ALB$A$equation2
#>                mean         sd 5% quantile 95% quantile
#> lag1_var1 -3.489234  0.6766011   -4.354866    -2.806352
#> lag1_var2  4.858919  0.5197574    4.278828     5.406216
#> lag1_var3  7.094315  0.7810578    6.131457     7.882793
#> lag1_var4 -5.796159  0.7752213   -6.602124    -4.860845
#> const      8.418724 14.1864859   -6.446161    26.157190
#> 
#> $ALB$A$equation3
#>                mean        sd 5% quantile 95% quantile
#> lag1_var1  2.532284  1.063908    1.210126     3.528763
#> lag1_var2 -3.062253  1.405690   -4.070698    -1.137900
#> lag1_var3 -4.390450  2.476307   -6.234374    -1.005255
#> lag1_var4  4.311895  2.059465    1.494075     5.870500
#> const      3.639306 11.668056   -4.262664    19.541851
#> 
#> $ALB$A$equation4
#>                   mean        sd 5% quantile 95% quantile
#> lag1_var1  0.213500502  1.217241   -1.347249     1.428493
#> lag1_var2 -0.080526422  1.336047   -1.002520     1.750552
#> lag1_var3 -0.004718148  2.499604   -1.669731     3.429913
#> lag1_var4  0.875370725  2.012951   -1.890817     2.266866
#> const      3.814572101 10.413278   -6.661226    14.973852
#> 
#> 
#> 
#> $ARE
#> $ARE$Sigma
#> $ARE$Sigma$equation1
#>                  mean         sd 5% quantile 95% quantile
#> Sigma[1,1] 0.01015687 0.00297673  0.00628476   0.01297308
#> 
#> $ARE$Sigma$equation2
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[2,1] -0.01328703 0.00897893 -0.02272233   -0.0021893
#> Sigma[2,2]  0.14230407 0.02568791  0.11570276    0.1746005
#> 
#> $ARE$Sigma$equation3
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[3,1]  0.03795118 0.01696096  0.02369078   0.06069523
#> Sigma[3,2] -0.18031514 0.05615531 -0.25202417  -0.12605440
#> Sigma[3,3]  0.45511713 0.14670603  0.32006545   0.64682452
#> 
#> $ARE$Sigma$equation4
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[4,1]  0.02891644 0.01183231  0.01792412   0.04434126
#> Sigma[4,2] -0.07018043 0.03674233 -0.11808156  -0.03700143
#> Sigma[4,3]  0.31720568 0.10420018  0.21108902   0.44941363
#> Sigma[4,4]  0.26711909 0.08133842  0.17913496   0.36398550
#> 
#> 
#> $ARE$A
#> $ARE$A$equation1
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1  0.9852596 0.04747217   0.9204773     1.019772
#> lag1_var2 -0.8857620 0.07923729  -0.9617789    -0.793329
#> lag1_var3 -1.1627453 0.06514790  -1.2481168    -1.103511
#> lag1_var4  1.1419182 0.06999885   1.0842840     1.236646
#> const      2.0078567 0.81683538   1.0048374     2.884372
#> 
#> $ARE$A$equation2
#>                 mean        sd 5% quantile 95% quantile
#> lag1_var1  0.2015970 0.2033726 -0.03378979   0.42995657
#> lag1_var2  0.3625828 0.2046322  0.12865852   0.60227812
#> lag1_var3 -0.2358660 0.2556469 -0.54346892   0.02892482
#> lag1_var4  0.2326629 0.2637937 -0.03494080   0.53144016
#> const     -3.9443736 2.9346243 -6.35377646  -0.19049027
#> 
#> $ARE$A$equation3
#>                 mean        sd 5% quantile 95% quantile
#> lag1_var1  0.2721741 0.6417852  -0.5080471    0.8787956
#> lag1_var2 -1.5224060 0.5907869  -2.2981170   -0.9936544
#> lag1_var3 -1.1742927 0.6679876  -2.0899654   -0.6997962
#> lag1_var4  1.9468618 0.6585135   1.3797941    2.8153253
#> const     10.0726796 9.1645745   0.4438359   20.6252563
#> 
#> $ARE$A$equation4
#>                mean        sd 5% quantile 95% quantile
#> lag1_var1  0.420213 0.7177969  -0.3516303    1.1562324
#> lag1_var2 -1.274005 0.5420638  -1.9506210   -0.7151177
#> lag1_var3 -1.390023 0.5877936  -2.1846052   -0.9594907
#> lag1_var4  2.195038 0.5267801   1.7948486    2.9150305
#> const      4.556889 9.3102913  -5.5426922   14.0866798
#> 
#> 
#> 
#> $ARG
#> $ARG$Sigma
#> $ARG$Sigma$equation1
#>                   mean          sd 5% quantile 95% quantile
#> Sigma[1,1] 0.005103857 0.002119963 0.003216289  0.007910083
#> 
#> $ARG$Sigma$equation2
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[2,1] -0.08128757 0.03680575  -0.1258673  -0.04213103
#> Sigma[2,2]  3.41520035 1.15371222   2.3661780   4.89939987
#> 
#> $ARG$Sigma$equation3
#>                   mean         sd 5% quantile 95% quantile
#> Sigma[3,1]  0.07413511 0.02034973  0.05403004    0.0963121
#> Sigma[3,2] -2.58462031 0.71817152 -3.46722460   -1.8752757
#> Sigma[3,3]  2.71239216 0.34287006  2.35855289    3.0725640
#> 
#> $ARG$Sigma$equation4
#>                  mean          sd 5% quantile 95% quantile
#> Sigma[4,1]  0.0287917 0.008033342  0.01975514   0.03818099
#> Sigma[4,2] -0.5952482 0.147768044 -0.74628427  -0.43750428
#> Sigma[4,3]  1.3079351 0.198830491  1.09226862   1.50127367
#> Sigma[4,4]  1.0794914 0.209278046  0.86997543   1.34535666
#> 
#> 
#> $ARG$A
#> $ARG$A$equation1
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1  0.9391390 0.06474795   0.8720649    1.0167766
#> lag1_var2  0.2717269 0.07850607   0.1829268    0.3617183
#> lag1_var3  0.4422313 0.12953779   0.2965005    0.5929919
#> lag1_var4 -0.4111026 0.11780038  -0.5464193   -0.2764494
#> const     -0.2539398 1.20774595  -1.7508025    1.0440643
#> 
#> $ARG$A$equation2
#>                  mean         sd 5% quantile 95% quantile
#> lag1_var1 -1.37823371  1.0714089  -2.5943901   -0.3756158
#> lag1_var2  0.08952226  0.7850127  -0.6504304    0.9663270
#> lag1_var3 -0.94351448  1.4104540  -2.2947439    0.6140563
#> lag1_var4  1.44424176  0.9993170   0.3004707    2.4076186
#> const     10.10621312 10.7499606  -2.4820113   20.9909336
#> 
#> $ARG$A$equation3
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1  0.7001508  0.4939051   0.1188289     1.248799
#> lag1_var2  2.8777674  0.6633487   1.9899993     3.417486
#> lag1_var3  5.3106551  1.1042644   3.8225342     6.172427
#> lag1_var4 -4.7217657  1.1490407  -5.7241037    -3.209022
#> const      3.6425928 12.1120860 -10.4729248    15.712324
#> 
#> $ARG$A$equation4
#>                 mean         sd 5% quantile 95% quantile
#> lag1_var1 -0.1607312  0.7043634  -0.9494391    0.6232343
#> lag1_var2  3.2256490  1.1024800   1.9349234    4.3785696
#> lag1_var3  5.2025452  1.8468812   3.0294285    7.1278721
#> lag1_var4 -4.2205533  1.7224182  -6.0864813   -2.2268985
#> const      4.9352412 15.5461656  -7.3121923   24.6392485
#> 
#> 
#> 
```
