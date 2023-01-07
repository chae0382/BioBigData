#Countgap 함수를 1000번 호출하여 계산한 평균 gap : 43.486
R version 4.2.1 (2022-06-23) -- "Funny-Looking Kid"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

During startup - Warning messages:
1: Setting LC_CTYPE failed, using "C" 
2: Setting LC_COLLATE failed, using "C" 
3: Setting LC_TIME failed, using "C" 
4: Setting LC_MESSAGES failed, using "C" 
5: Setting LC_MONETARY failed, using "C" 
[R.app GUI 1.79 (8095) x86_64-apple-darwin17.0]

WARNING: You're using a non-UTF8 locale, therefore only ASCII characters will work.
Please read R for Mac OS X FAQ (see Help) section 9 and adjust your system preferences accordingly.
[Workspace restored from /Users/chaeyeon/.RData]
[History restored from /Users/chaeyeon/.Rapp.history]

2022-10-08 22:57:56.911 R[72490:874590] TSM AdjustCapsLockLEDForKeyTransitionHandling - _ISSetPhysicalKeyboardCapsLockLED Inhibit
> #gap의 개수를 계산하는 함수 정의
> Countgap=function(N, l, m)
+ {
+ X = sample(1:N, size=m, replace=T)
+ sort_X = sort(X)
+ gap = 0
+ num = sum(sort_X< (N-l+2))
+ for(i in 1:(num-1)){
+ if(sort_X[i]+l-1 < sort_X[i+1]){
+ gap = gap + 1
+ }
+ }
+ return(gap)
+ }
> #위의 함수를 1000번 호출하여 평균 gap 구하기
> sum_of_gaps = 0
> for(i in 1:1000){
+ gap = Countgap(10000,10,5000)
+ sum_of_gaps = sum_of_gaps + gap
+ }
> average_number_of_gaps = sum_of_gaps/1000
> average_number_of_gaps
[1] 43.486
> 