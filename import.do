
clear

fdause "~/Desktop/tmpah2/intercultural/allwave1.xpt" // data from wave 1 inhome + inschool surveys

/// keep relevant variables only
keep aid scid commid bio_sex h1gi1m h1gi1y h1gi3 h1gi4 h1gi6a h1gi6b h1gi6c h1gi6d h1gi6e h1gi20 h1rm1 h1rm4 h1rm5 h1rf1 h1rf4 h1rf5 h1wp1 h1wp2 h1wp3 h1wp4 h1wp5 h1wp6 h1wp7 h1wp8 h1wp9 h1wp10 h1wp13 h1wp14 h1wp17h h1wp17i h1wp17j h1wp18h h1wp18i h1wp18j fr_flag h1pr3 h1pr4 h1nb1 h1nb2 h1nb3 h1nb4 h1nb5 h1nb6 h1nb7 h1re1 sqid sschlcde s3 s4 s6a s6b s6c s6d s6e s9 s11 s12 s14 s15 s16 s17 s18 s19 s20 s21 s22 s44a1-s44 pa12 pa13 pb9 h1da4 h1da5 h1da6 h1da7

// sort by unique student id
sort aid

save ~/Desktop/tmpah2/intercultural/newdata, replace

fdause "~/Desktop/tmpah2/intercultural/hfriend1.xpt", clear // data on friendship relations, inhome survey, wave1
sort aid
merge aid using ~/Desktop/tmpah2/intercultural/newdata

keep if _merge==3
drop _merge
sort aid

save ~/Desktop/tmpah2/intercultural/newdata, replace

fdause "~/Desktop/tmpah2/intercultural/Spatial.xpt", clear // data on normalized geographical locations
sort aid

merge aid using ~/Desktop/tmpah2/intercultural/newdata

keep if _merge==3
drop _merge
sort aid

destring _all, replace force

keep if scid==1 | scid==2 | scid==3 | scid==7 | scid==8 | scid==28 | scid==58 | scid==77 | scid==81 | scid==88 | scid==106 | scid==115 | scid==126 | scid==175 | scid==194 | scid==369   //keep saturated schools

mvencode _all, mv(-1) // encore missing values

save ~/Desktop/tmpah2/intercultural/newdata, replace
