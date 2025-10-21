# BayesCOOP - Bayesian Cooperative Learning for Multimodal Integration

The repository houses the **`BayesCOOP`** R package for multimodal
integration for prediction of continuous outcomes.

## Dependencies

`BayesCOOP` requires the following `R` package: `devtools` (for
installation only). Please install it before installing `BayesCOOP`,
which can be done as follows (execute from within a fresh R session):

    install.packages("devtools")
    library(devtools)

## Installation

Once the dependencies are installed, `BayesCOOP` can be loaded using the
following command:

    devtools::install_github("himelmallick/BayesCOOP", quiet = TRUE)
    library(BayesCOOP)

## Example Implementation

### Loading the StelzerDOS real dataset

    data_train = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS.RData"))); rm(pcl)
    data_test = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS_valid.RData"))); rm(pcl)

### Pre-processing the longitudinal data by considering only baseline observations

#### Remove metabolomics from the train set to match with validation

    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr", repos = "https://cloud.r-project.org")
    }
    library(dplyr)

    data_train$feature_metadata = data_train$feature_metadata %>% dplyr::filter(featureType!='Metabolomics')
    data_train$feature_table = data_train$feature_table[rownames(data_train$feature_metadata),]

#### Consider only baseline observations for the train set

    positions = grep("A", colnames(data_train$feature_table), ignore.case = TRUE)
    data_train$feature_table = data_train$feature_table[, positions]
    data_train$sample_metadata = data_train$sample_metadata[positions, ]
    rm(positions)

#### Consider only baseline observations for the validation set

    positions = grep("G1", colnames(data_test$feature_table))
    data_test$feature_table = data_test$feature_table[, positions]
    data_test$sample_metadata = data_test$sample_metadata[positions, ]
    rm(positions)

### Implementing BayesCOOP

    set.seed(1)
    result = BayesCOOP::bayesCoop(data_train, data_test, family = "gaussian", 
                          ss = c(0.05, 1), group = TRUE,
                          bb = TRUE, alpha_dirich = 1, 
                          bbiters = 1100, bbburn = 100, maxit = 100,
                          filter = TRUE, abd_thresh = 0, prev_thresh = 0.1,
                          Warning = TRUE, verbose = TRUE, control = list())
    ## Implementing BayesCoop...
    ## iter = 1
    ## iter = 2
    ## iter = 3
    ## iter = 4
    ## iter = 5
    ## iter = 6
    ## iter = 7
    ## iter = 8
    ## iter = 9
    ## iter = 10
    ## iter = 11
    ## iter = 12
    ## iter = 13
    ## iter = 14
    ## iter = 15
    ## iter = 16
    ## iter = 17
    ## iter = 18
    ## iter = 19
    ## iter = 20
    ## iter = 21
    ## iter = 22
    ## iter = 23
    ## iter = 24
    ## iter = 25
    ## iter = 26
    ## iter = 27
    ## iter = 28
    ## iter = 29
    ## iter = 30
    ## iter = 31
    ## iter = 32
    ## iter = 33
    ## iter = 34
    ## iter = 35
    ## iter = 36
    ## iter = 37
    ## iter = 38
    ## iter = 39
    ## iter = 40
    ## iter = 41
    ## iter = 42
    ## iter = 43
    ## iter = 44
    ## iter = 45
    ## iter = 46
    ## iter = 47
    ## iter = 48
    ## iter = 49
    ## iter = 50
    ## iter = 51
    ## iter = 52
    ## iter = 53
    ## iter = 54
    ## iter = 55
    ## iter = 56
    ## iter = 57
    ## iter = 58
    ## iter = 59
    ## iter = 60
    ## iter = 61
    ## iter = 62
    ## iter = 63
    ## iter = 64
    ## iter = 65
    ## iter = 66
    ## iter = 67
    ## iter = 68
    ## iter = 69
    ## iter = 70
    ## iter = 71
    ## iter = 72
    ## iter = 73
    ## iter = 74
    ## iter = 75
    ## iter = 76
    ## iter = 77
    ## iter = 78
    ## iter = 79
    ## iter = 80
    ## iter = 81
    ## iter = 82
    ## iter = 83
    ## iter = 84
    ## iter = 85
    ## iter = 86
    ## iter = 87
    ## iter = 88
    ## iter = 89
    ## iter = 90
    ## iter = 91
    ## iter = 92
    ## iter = 93
    ## iter = 94
    ## iter = 95
    ## iter = 96
    ## iter = 97
    ## iter = 98
    ## iter = 99
    ## iter = 100
    ## iter = 101
    ## iter = 102
    ## iter = 103
    ## iter = 104
    ## iter = 105
    ## iter = 106
    ## iter = 107
    ## iter = 108
    ## iter = 109
    ## iter = 110
    ## iter = 111
    ## iter = 112
    ## iter = 113
    ## iter = 114
    ## iter = 115
    ## iter = 116
    ## iter = 117
    ## iter = 118
    ## iter = 119
    ## iter = 120
    ## iter = 121
    ## iter = 122
    ## iter = 123
    ## iter = 124
    ## iter = 125
    ## iter = 126
    ## iter = 127
    ## iter = 128
    ## iter = 129
    ## iter = 130
    ## iter = 131
    ## iter = 132
    ## iter = 133
    ## iter = 134
    ## iter = 135
    ## iter = 136
    ## iter = 137
    ## iter = 138
    ## iter = 139
    ## iter = 140
    ## iter = 141
    ## iter = 142
    ## iter = 143
    ## iter = 144
    ## iter = 145
    ## iter = 146
    ## iter = 147
    ## iter = 148
    ## iter = 149
    ## iter = 150
    ## iter = 151
    ## iter = 152
    ## iter = 153
    ## iter = 154
    ## iter = 155
    ## iter = 156
    ## iter = 157
    ## iter = 158
    ## iter = 159
    ## iter = 160
    ## iter = 161
    ## iter = 162
    ## iter = 163
    ## iter = 164
    ## iter = 165
    ## iter = 166
    ## iter = 167
    ## iter = 168
    ## iter = 169
    ## iter = 170
    ## iter = 171
    ## iter = 172
    ## iter = 173
    ## iter = 174
    ## iter = 175
    ## iter = 176
    ## iter = 177
    ## iter = 178
    ## iter = 179
    ## iter = 180
    ## iter = 181
    ## iter = 182
    ## iter = 183
    ## iter = 184
    ## iter = 185
    ## iter = 186
    ## iter = 187
    ## iter = 188
    ## iter = 189
    ## iter = 190
    ## iter = 191
    ## iter = 192
    ## iter = 193
    ## iter = 194
    ## iter = 195
    ## iter = 196
    ## iter = 197
    ## iter = 198
    ## iter = 199
    ## iter = 200
    ## iter = 201
    ## iter = 202
    ## iter = 203
    ## iter = 204
    ## iter = 205
    ## iter = 206
    ## iter = 207
    ## iter = 208
    ## iter = 209
    ## iter = 210
    ## iter = 211
    ## iter = 212
    ## iter = 213
    ## iter = 214
    ## iter = 215
    ## iter = 216
    ## iter = 217
    ## iter = 218
    ## iter = 219
    ## iter = 220
    ## iter = 221
    ## iter = 222
    ## iter = 223
    ## iter = 224
    ## iter = 225
    ## iter = 226
    ## iter = 227
    ## iter = 228
    ## iter = 229
    ## iter = 230
    ## iter = 231
    ## iter = 232
    ## iter = 233
    ## iter = 234
    ## iter = 235
    ## iter = 236
    ## iter = 237
    ## iter = 238
    ## iter = 239
    ## iter = 240
    ## iter = 241
    ## iter = 242
    ## iter = 243
    ## iter = 244
    ## iter = 245
    ## iter = 246
    ## iter = 247
    ## iter = 248
    ## iter = 249
    ## iter = 250
    ## iter = 251
    ## iter = 252
    ## iter = 253
    ## iter = 254
    ## iter = 255
    ## iter = 256
    ## iter = 257
    ## iter = 258
    ## iter = 259
    ## iter = 260
    ## iter = 261
    ## iter = 262
    ## iter = 263
    ## iter = 264
    ## iter = 265
    ## iter = 266
    ## iter = 267
    ## iter = 268
    ## iter = 269
    ## iter = 270
    ## iter = 271
    ## iter = 272
    ## iter = 273
    ## iter = 274
    ## iter = 275
    ## iter = 276
    ## iter = 277
    ## iter = 278
    ## iter = 279
    ## iter = 280
    ## iter = 281
    ## iter = 282
    ## iter = 283
    ## iter = 284
    ## iter = 285
    ## iter = 286
    ## iter = 287
    ## iter = 288
    ## iter = 289
    ## iter = 290
    ## iter = 291
    ## iter = 292
    ## iter = 293
    ## iter = 294
    ## iter = 295
    ## iter = 296
    ## iter = 297
    ## iter = 298
    ## iter = 299
    ## iter = 300
    ## iter = 301
    ## iter = 302
    ## iter = 303
    ## iter = 304
    ## iter = 305
    ## iter = 306
    ## iter = 307
    ## iter = 308
    ## iter = 309
    ## iter = 310
    ## iter = 311
    ## iter = 312
    ## iter = 313
    ## iter = 314
    ## iter = 315
    ## iter = 316
    ## iter = 317
    ## iter = 318
    ## iter = 319
    ## iter = 320
    ## iter = 321
    ## iter = 322
    ## iter = 323
    ## iter = 324
    ## iter = 325
    ## iter = 326
    ## iter = 327
    ## iter = 328
    ## iter = 329
    ## iter = 330
    ## iter = 331
    ## iter = 332
    ## iter = 333
    ## iter = 334
    ## iter = 335
    ## iter = 336
    ## iter = 337
    ## iter = 338
    ## iter = 339
    ## iter = 340
    ## iter = 341
    ## iter = 342
    ## iter = 343
    ## iter = 344
    ## iter = 345
    ## iter = 346
    ## iter = 347
    ## iter = 348
    ## iter = 349
    ## iter = 350
    ## iter = 351
    ## iter = 352
    ## iter = 353
    ## iter = 354
    ## iter = 355
    ## iter = 356
    ## iter = 357
    ## iter = 358
    ## iter = 359
    ## iter = 360
    ## iter = 361
    ## iter = 362
    ## iter = 363
    ## iter = 364
    ## iter = 365
    ## iter = 366
    ## iter = 367
    ## iter = 368
    ## iter = 369
    ## iter = 370
    ## iter = 371
    ## iter = 372
    ## iter = 373
    ## iter = 374
    ## iter = 375
    ## iter = 376
    ## iter = 377
    ## iter = 378
    ## iter = 379
    ## iter = 380
    ## iter = 381
    ## iter = 382
    ## iter = 383
    ## iter = 384
    ## iter = 385
    ## iter = 386
    ## iter = 387
    ## iter = 388
    ## iter = 389
    ## iter = 390
    ## iter = 391
    ## iter = 392
    ## iter = 393
    ## iter = 394
    ## iter = 395
    ## iter = 396
    ## iter = 397
    ## iter = 398
    ## iter = 399
    ## iter = 400
    ## iter = 401
    ## iter = 402
    ## iter = 403
    ## iter = 404
    ## iter = 405
    ## iter = 406
    ## iter = 407
    ## iter = 408
    ## iter = 409
    ## iter = 410
    ## iter = 411
    ## iter = 412
    ## iter = 413
    ## iter = 414
    ## iter = 415
    ## iter = 416
    ## iter = 417
    ## iter = 418
    ## iter = 419
    ## iter = 420
    ## iter = 421
    ## iter = 422
    ## iter = 423
    ## iter = 424
    ## iter = 425
    ## iter = 426
    ## iter = 427
    ## iter = 428
    ## iter = 429
    ## iter = 430
    ## iter = 431
    ## iter = 432
    ## iter = 433
    ## iter = 434
    ## iter = 435
    ## iter = 436
    ## iter = 437
    ## iter = 438
    ## iter = 439
    ## iter = 440
    ## iter = 441
    ## iter = 442
    ## iter = 443
    ## iter = 444
    ## iter = 445
    ## iter = 446
    ## iter = 447
    ## iter = 448
    ## iter = 449
    ## iter = 450
    ## iter = 451
    ## iter = 452
    ## iter = 453
    ## iter = 454
    ## iter = 455
    ## iter = 456
    ## iter = 457
    ## iter = 458
    ## iter = 459
    ## iter = 460
    ## iter = 461
    ## iter = 462
    ## iter = 463
    ## iter = 464
    ## iter = 465
    ## iter = 466
    ## iter = 467
    ## iter = 468
    ## iter = 469
    ## iter = 470
    ## iter = 471
    ## iter = 472
    ## iter = 473
    ## iter = 474
    ## iter = 475
    ## iter = 476
    ## iter = 477
    ## iter = 478
    ## iter = 479
    ## iter = 480
    ## iter = 481
    ## iter = 482
    ## iter = 483
    ## iter = 484
    ## iter = 485
    ## iter = 486
    ## iter = 487
    ## iter = 488
    ## iter = 489
    ## iter = 490
    ## iter = 491
    ## iter = 492
    ## iter = 493
    ## iter = 494
    ## iter = 495
    ## iter = 496
    ## iter = 497
    ## iter = 498
    ## iter = 499
    ## iter = 500
    ## iter = 501
    ## iter = 502
    ## iter = 503
    ## iter = 504
    ## iter = 505
    ## iter = 506
    ## iter = 507
    ## iter = 508
    ## iter = 509
    ## iter = 510
    ## iter = 511
    ## iter = 512
    ## iter = 513
    ## iter = 514
    ## iter = 515
    ## iter = 516
    ## iter = 517
    ## iter = 518
    ## iter = 519
    ## iter = 520
    ## iter = 521
    ## iter = 522
    ## iter = 523
    ## iter = 524
    ## iter = 525
    ## iter = 526
    ## iter = 527
    ## iter = 528
    ## iter = 529
    ## iter = 530
    ## iter = 531
    ## iter = 532
    ## iter = 533
    ## iter = 534
    ## iter = 535
    ## iter = 536
    ## iter = 537
    ## iter = 538
    ## iter = 539
    ## iter = 540
    ## iter = 541
    ## iter = 542
    ## iter = 543
    ## iter = 544
    ## iter = 545
    ## iter = 546
    ## iter = 547
    ## iter = 548
    ## iter = 549
    ## iter = 550
    ## iter = 551
    ## iter = 552
    ## iter = 553
    ## iter = 554
    ## iter = 555
    ## iter = 556
    ## iter = 557
    ## iter = 558
    ## iter = 559
    ## iter = 560
    ## iter = 561
    ## iter = 562
    ## iter = 563
    ## iter = 564
    ## iter = 565
    ## iter = 566
    ## iter = 567
    ## iter = 568
    ## iter = 569
    ## iter = 570
    ## iter = 571
    ## iter = 572
    ## iter = 573
    ## iter = 574
    ## iter = 575
    ## iter = 576
    ## iter = 577
    ## iter = 578
    ## iter = 579
    ## iter = 580
    ## iter = 581
    ## iter = 582
    ## iter = 583
    ## iter = 584
    ## iter = 585
    ## iter = 586
    ## iter = 587
    ## iter = 588
    ## iter = 589
    ## iter = 590
    ## iter = 591
    ## iter = 592
    ## iter = 593
    ## iter = 594
    ## iter = 595
    ## iter = 596
    ## iter = 597
    ## iter = 598
    ## iter = 599
    ## iter = 600
    ## iter = 601
    ## iter = 602
    ## iter = 603
    ## iter = 604
    ## iter = 605
    ## iter = 606
    ## iter = 607
    ## iter = 608
    ## iter = 609
    ## iter = 610
    ## iter = 611
    ## iter = 612
    ## iter = 613
    ## iter = 614
    ## iter = 615
    ## iter = 616
    ## iter = 617
    ## iter = 618
    ## iter = 619
    ## iter = 620
    ## iter = 621
    ## iter = 622
    ## iter = 623
    ## iter = 624
    ## iter = 625
    ## iter = 626
    ## iter = 627
    ## iter = 628
    ## iter = 629
    ## iter = 630
    ## iter = 631
    ## iter = 632
    ## iter = 633
    ## iter = 634
    ## iter = 635
    ## iter = 636
    ## iter = 637
    ## iter = 638
    ## iter = 639
    ## iter = 640
    ## iter = 641
    ## iter = 642
    ## iter = 643
    ## iter = 644
    ## iter = 645
    ## iter = 646
    ## iter = 647
    ## iter = 648
    ## iter = 649
    ## iter = 650
    ## iter = 651
    ## iter = 652
    ## iter = 653
    ## iter = 654
    ## iter = 655
    ## iter = 656
    ## iter = 657
    ## iter = 658
    ## iter = 659
    ## iter = 660
    ## iter = 661
    ## iter = 662
    ## iter = 663
    ## iter = 664
    ## iter = 665
    ## iter = 666
    ## iter = 667
    ## iter = 668
    ## iter = 669
    ## iter = 670
    ## iter = 671
    ## iter = 672
    ## iter = 673
    ## iter = 674
    ## iter = 675
    ## iter = 676
    ## iter = 677
    ## iter = 678
    ## iter = 679
    ## iter = 680
    ## iter = 681
    ## iter = 682
    ## iter = 683
    ## iter = 684
    ## iter = 685
    ## iter = 686
    ## iter = 687
    ## iter = 688
    ## iter = 689
    ## iter = 690
    ## iter = 691
    ## iter = 692
    ## iter = 693
    ## iter = 694
    ## iter = 695
    ## iter = 696
    ## iter = 697
    ## iter = 698
    ## iter = 699
    ## iter = 700
    ## iter = 701
    ## iter = 702
    ## iter = 703
    ## iter = 704
    ## iter = 705
    ## iter = 706
    ## iter = 707
    ## iter = 708
    ## iter = 709
    ## iter = 710
    ## iter = 711
    ## iter = 712
    ## iter = 713
    ## iter = 714
    ## iter = 715
    ## iter = 716
    ## iter = 717
    ## iter = 718
    ## iter = 719
    ## iter = 720
    ## iter = 721
    ## iter = 722
    ## iter = 723
    ## iter = 724
    ## iter = 725
    ## iter = 726
    ## iter = 727
    ## iter = 728
    ## iter = 729
    ## iter = 730
    ## iter = 731
    ## iter = 732
    ## iter = 733
    ## iter = 734
    ## iter = 735
    ## iter = 736
    ## iter = 737
    ## iter = 738
    ## iter = 739
    ## iter = 740
    ## iter = 741
    ## iter = 742
    ## iter = 743
    ## iter = 744
    ## iter = 745
    ## iter = 746
    ## iter = 747
    ## iter = 748
    ## iter = 749
    ## iter = 750
    ## iter = 751
    ## iter = 752
    ## iter = 753
    ## iter = 754
    ## iter = 755
    ## iter = 756
    ## iter = 757
    ## iter = 758
    ## iter = 759
    ## iter = 760
    ## iter = 761
    ## iter = 762
    ## iter = 763
    ## iter = 764
    ## iter = 765
    ## iter = 766
    ## iter = 767
    ## iter = 768
    ## iter = 769
    ## iter = 770
    ## iter = 771
    ## iter = 772
    ## iter = 773
    ## iter = 774
    ## iter = 775
    ## iter = 776
    ## iter = 777
    ## iter = 778
    ## iter = 779
    ## iter = 780
    ## iter = 781
    ## iter = 782
    ## iter = 783
    ## iter = 784
    ## iter = 785
    ## iter = 786
    ## iter = 787
    ## iter = 788
    ## iter = 789
    ## iter = 790
    ## iter = 791
    ## iter = 792
    ## iter = 793
    ## iter = 794
    ## iter = 795
    ## iter = 796
    ## iter = 797
    ## iter = 798
    ## iter = 799
    ## iter = 800
    ## iter = 801
    ## iter = 802
    ## iter = 803
    ## iter = 804
    ## iter = 805
    ## iter = 806
    ## iter = 807
    ## iter = 808
    ## iter = 809
    ## iter = 810
    ## iter = 811
    ## iter = 812
    ## iter = 813
    ## iter = 814
    ## iter = 815
    ## iter = 816
    ## iter = 817
    ## iter = 818
    ## iter = 819
    ## iter = 820
    ## iter = 821
    ## iter = 822
    ## iter = 823
    ## iter = 824
    ## iter = 825
    ## iter = 826
    ## iter = 827
    ## iter = 828
    ## iter = 829
    ## iter = 830
    ## iter = 831
    ## iter = 832
    ## iter = 833
    ## iter = 834
    ## iter = 835
    ## iter = 836
    ## iter = 837
    ## iter = 838
    ## iter = 839
    ## iter = 840
    ## iter = 841
    ## iter = 842
    ## iter = 843
    ## iter = 844
    ## iter = 845
    ## iter = 846
    ## iter = 847
    ## iter = 848
    ## iter = 849
    ## iter = 850
    ## iter = 851
    ## iter = 852
    ## iter = 853
    ## iter = 854
    ## iter = 855
    ## iter = 856
    ## iter = 857
    ## iter = 858
    ## iter = 859
    ## iter = 860
    ## iter = 861
    ## iter = 862
    ## iter = 863
    ## iter = 864
    ## iter = 865
    ## iter = 866
    ## iter = 867
    ## iter = 868
    ## iter = 869
    ## iter = 870
    ## iter = 871
    ## iter = 872
    ## iter = 873
    ## iter = 874
    ## iter = 875
    ## iter = 876
    ## iter = 877
    ## iter = 878
    ## iter = 879
    ## iter = 880
    ## iter = 881
    ## iter = 882
    ## iter = 883
    ## iter = 884
    ## iter = 885
    ## iter = 886
    ## iter = 887
    ## iter = 888
    ## iter = 889
    ## iter = 890
    ## iter = 891
    ## iter = 892
    ## iter = 893
    ## iter = 894
    ## iter = 895
    ## iter = 896
    ## iter = 897
    ## iter = 898
    ## iter = 899
    ## iter = 900
    ## iter = 901
    ## iter = 902
    ## iter = 903
    ## iter = 904
    ## iter = 905
    ## iter = 906
    ## iter = 907
    ## iter = 908
    ## iter = 909
    ## iter = 910
    ## iter = 911
    ## iter = 912
    ## iter = 913
    ## iter = 914
    ## iter = 915
    ## iter = 916
    ## iter = 917
    ## iter = 918
    ## iter = 919
    ## iter = 920
    ## iter = 921
    ## iter = 922
    ## iter = 923
    ## iter = 924
    ## iter = 925
    ## iter = 926
    ## iter = 927
    ## iter = 928
    ## iter = 929
    ## iter = 930
    ## iter = 931
    ## iter = 932
    ## iter = 933
    ## iter = 934
    ## iter = 935
    ## iter = 936
    ## iter = 937
    ## iter = 938
    ## iter = 939
    ## iter = 940
    ## iter = 941
    ## iter = 942
    ## iter = 943
    ## iter = 944
    ## iter = 945
    ## iter = 946
    ## iter = 947
    ## iter = 948
    ## iter = 949
    ## iter = 950
    ## iter = 951
    ## iter = 952
    ## iter = 953
    ## iter = 954
    ## iter = 955
    ## iter = 956
    ## iter = 957
    ## iter = 958
    ## iter = 959
    ## iter = 960
    ## iter = 961
    ## iter = 962
    ## iter = 963
    ## iter = 964
    ## iter = 965
    ## iter = 966
    ## iter = 967
    ## iter = 968
    ## iter = 969
    ## iter = 970
    ## iter = 971
    ## iter = 972
    ## iter = 973
    ## iter = 974
    ## iter = 975
    ## iter = 976
    ## iter = 977
    ## iter = 978
    ## iter = 979
    ## iter = 980
    ## iter = 981
    ## iter = 982
    ## iter = 983
    ## iter = 984
    ## iter = 985
    ## iter = 986
    ## iter = 987
    ## iter = 988
    ## iter = 989
    ## iter = 990
    ## iter = 991
    ## iter = 992
    ## iter = 993
    ## iter = 994
    ## iter = 995
    ## iter = 996
    ## iter = 997
    ## iter = 998
    ## iter = 999
    ## iter = 1000
    ## iter = 1001
    ## iter = 1002
    ## iter = 1003
    ## iter = 1004
    ## iter = 1005
    ## iter = 1006
    ## iter = 1007
    ## iter = 1008
    ## iter = 1009
    ## iter = 1010
    ## iter = 1011
    ## iter = 1012
    ## iter = 1013
    ## iter = 1014
    ## iter = 1015
    ## iter = 1016
    ## iter = 1017
    ## iter = 1018
    ## iter = 1019
    ## iter = 1020
    ## iter = 1021
    ## iter = 1022
    ## iter = 1023
    ## iter = 1024
    ## iter = 1025
    ## iter = 1026
    ## iter = 1027
    ## iter = 1028
    ## iter = 1029
    ## iter = 1030
    ## iter = 1031
    ## iter = 1032
    ## iter = 1033
    ## iter = 1034
    ## iter = 1035
    ## iter = 1036
    ## iter = 1037
    ## iter = 1038
    ## iter = 1039
    ## iter = 1040
    ## iter = 1041
    ## iter = 1042
    ## iter = 1043
    ## iter = 1044
    ## iter = 1045
    ## iter = 1046
    ## iter = 1047
    ## iter = 1048
    ## iter = 1049
    ## iter = 1050
    ## iter = 1051
    ## iter = 1052
    ## iter = 1053
    ## iter = 1054
    ## iter = 1055
    ## iter = 1056
    ## iter = 1057
    ## iter = 1058
    ## iter = 1059
    ## iter = 1060
    ## iter = 1061
    ## iter = 1062
    ## iter = 1063
    ## iter = 1064
    ## iter = 1065
    ## iter = 1066
    ## iter = 1067
    ## iter = 1068
    ## iter = 1069
    ## iter = 1070
    ## iter = 1071
    ## iter = 1072
    ## iter = 1073
    ## iter = 1074
    ## iter = 1075
    ## iter = 1076
    ## iter = 1077
    ## iter = 1078
    ## iter = 1079
    ## iter = 1080
    ## iter = 1081
    ## iter = 1082
    ## iter = 1083
    ## iter = 1084
    ## iter = 1085
    ## iter = 1086
    ## iter = 1087
    ## iter = 1088
    ## iter = 1089
    ## iter = 1090
    ## iter = 1091
    ## iter = 1092
    ## iter = 1093
    ## iter = 1094
    ## iter = 1095
    ## iter = 1096
    ## iter = 1097
    ## iter = 1098
    ## iter = 1099
    ## iter = 1100

### Implementing Cooperative Learning \[1\]

    if (!requireNamespace("multiview", quietly = TRUE)) {
      install.packages("multiview", repos = "https://cloud.r-project.org")
    }
    library(multiview)

    y_train <- BayesCOOP:::gen_datalist(data_train)$y
    y_test <- BayesCOOP:::gen_datalist(data_test)$y
    xList_train <- BayesCOOP:::gen_datalist(data_train)$xList
    xList_test <- BayesCOOP:::gen_datalist(data_test)$xList

    ## Filtering features
    xList_train <- lapply(xList_train, function(foo) BayesCOOP:::filter_features(foo, abd_thresh = 0, prev_thresh = 0.1))
    xList_test <- lapply(xList_test, function(foo) BayesCOOP:::filter_features(foo, abd_thresh = 0, prev_thresh = 0.1))

    ## Considering the common features between train and test set
    keep_features <- vector("list", length = length(xList_train))
    for(i in 1:length(keep_features)) {
        keep_features[[i]] <- intersect(names(xList_train[[i]]), names(xList_test[[i]]))
        xList_train[[i]] <- as.matrix(xList_train[[i]][, keep_features[[i]], drop = FALSE])
        xList_test[[i]] <- as.matrix(xList_test[[i]][, keep_features[[i]], drop = FALSE])
    }

    xList_train <- lapply(xList_train, function(foo) as.matrix(foo))
    xList_test <- lapply(xList_test, function(foo) as.matrix(foo)) ## need matrix form for multiview

    # Function to perform cross-validation for choosing alpha for cooperative learning application
    cvfit_alpha <- function(xlist, y, rho, alpha, fold) {
      cvfit <- cv.multiview(x_list = xlist, y = y, rho = rho,
                           alpha = alpha, nfolds = fold)
      lmin <- cvfit$lambda.min
      lmin_indx <- which(cvfit$lambda == lmin)
      cverr_min <- cvfit$cvm[lmin_indx]
      return(list(lmin = lmin, lmin_indx = lmin_indx, cverr_min = cverr_min,
                  alpha = alpha))
    }

    start.time <- Sys.time()
    cv_results_early <- cvfit_alpha(xlist = xList_train, y = y_train, rho = 0, alpha = 1, fold = 5)
    opt_lambda_early <- cv_results_early$lmin
    fit_multiview_early <- multiview(xList_train, y_train, lambda = opt_lambda_early, rho = 0, alpha = 1)
    stop.time <- Sys.time()
    time_multiview_early <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
    rm(start.time); rm(stop.time)
    beta_hat_early <- as.vector(fit_multiview_early$beta)
    ## prediction
    y_pred_early <- predict(fit_multiview_early, newx = xList_test, s = opt_lambda_early, rho = 0, alpha = 1)
    ## mean squared prediction error 
    mspe_multiview_early <- mean((y_pred_early - y_test) ^ 2)

    start.time <- Sys.time()
    cv_results_intermediate <- cvfit_alpha(xlist = xList_train, y = y_train, rho = 0.5, alpha = 1, fold = 5)
    opt_lambda_intermediate <- cv_results_intermediate$lmin
    fit_multiview_intermediate <- multiview(xList_train, y_train, lambda = opt_lambda_intermediate, rho = 0.5, alpha = 1)
    stop.time <- Sys.time()
    time_multiview_intermediate <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
    rm(start.time); rm(stop.time)
    beta_hat_intermediate <- as.vector(fit_multiview_intermediate$beta)
    ## prediction
    y_pred_intermediate <- predict(fit_multiview_intermediate, newx = xList_test, s = opt_lambda_intermediate, rho = 0.5, alpha = 1)
    ## mean squared prediction error 
    mspe_multiview_intermediate <- mean((y_pred_intermediate - y_test) ^ 2)

    start.time <- Sys.time()
    cv_results_late <- cvfit_alpha(xlist = xList_train, y = y_train, rho = 1, alpha = 1, fold = 5)
    opt_lambda_late <- cv_results_late$lmin
    fit_multiview_late <- multiview(xList_train, y_train, lambda = opt_lambda_late, rho = 1, alpha = 1)
    stop.time <- Sys.time()
    time_multiview_late <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
    rm(start.time); rm(stop.time)
    beta_hat_late <- as.vector(fit_multiview_late$beta)
    ## prediction
    y_pred_late <- predict(fit_multiview_late, newx = xList_test, s = opt_lambda_late, rho = 1, alpha = 1)
    ## mean squared prediction error 
    mspe_multiview_late <- mean((y_pred_late - y_test) ^ 2)

### Comparing BayesCOOP with Cooperative Learning

    ## BayesCOOP
    print(result$mspe)
    ## [1] 501.0142
    print(result$time)
    ## [1] 4.247

    ## Early Fusion
    print(mspe_multiview_early)
    ## [1] 869.5068
    print(time_multiview_early)
    ## [1] 0.003
    ## Intermediate Fusion
    print(mspe_multiview_intermediate)
    ## [1] 824.416
    print(time_multiview_intermediate)
    ## [1] 0.003
    ## Late Fusion
    print(mspe_multiview_late)
    ## [1] 948.7386
    print(time_multiview_late)
    ## [1] 0.005

## References

\[1\] Daisy Yi Ding, Shuangning Li, Balasubramanian Narasimhan, and
Robert Tibshirani (2022). “Cooperative Learning for Multiview Analysis.”
Proceedings of the National Academy of Sciences (USA), Vol. 119, No. 38,
e2202113119. DOI:
[10.1073/pnas.2202113119](https://www.pnas.org/doi/full/10.1073/pnas.2202113119)
