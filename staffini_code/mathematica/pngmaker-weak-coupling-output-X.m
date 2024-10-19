(* ::Package:: *)

(* ::Title::Closed:: *)
(*Plotter*)

filename=FileNames["*h5"][[1]]


(* find up to five characters after the var names, discard spurious. This means this works only for values of parameters up to 10000 *)
posFP=StringPosition[filename,"fP="] [[1]][[2]];
pos\[CapitalDelta]=StringPosition[filename,"Delta0="] [[1]][[2]];
fP=StringSplit[StringTake[filename,{posFP+1,posFP+6}],"_"][[1]]
Delta=StringSplit[StringTake[filename,{pos\[CapitalDelta]+1,pos\[CapitalDelta]+6}],"_"][[1]]

 posgrid = StringPosition[filename, "-grid_"] [[1]][[2]];
posrange = StringPosition[filename, "-range_"] [[1]][[2]];
fP = StringSplit[StringTake[filename, {posFP + 1, posFP + 6}], "_"][[
  1]]
Delta = StringSplit[
   StringTake[
    filename, {pos\[CapitalDelta] + 1, pos\[CapitalDelta] + 6}], 
   "_"][[1]]
grid = StringSplit[StringTake[filename, {posgrid - 9, posgrid - 6}], 
   "_"][[1]]
range = StringSplit[
   StringTake[filename, {posrange - 9, posrange - 7}], "_"][[1]]



t1 = Import[filename, {"Datasets", "/1/t"}];
r1 = Import[filename, {"Datasets", "/1/x"}];
imatom1 = Import[filename, {"Datasets", "/1/imatom"}];
reatom1 = Import[filename, {"Datasets", "/1/reatom"}];
modatomsq1 = Import[filename, {"Datasets", "/1/modatom"}];


declaredVariables={"t1", "r1", "imatom1", "reatom1", "imlight1", "relight1", "modatomsq1"}
tmax=Length[t1]
rmax=Length[r1]
{0,Last[r1]}

plotmin=1;
plotmax=tmax;
plotint=200;

imgsize=300;

(*pairGridPlot[plot1_,plot2_]:=
GraphicsGrid[Table[{ListPlot[Flatten[Table[{r1[[i]],plot1[[t,i]] } ,{i,1,rmax}],0], Joined->True, ImageSize->300,PlotRange->All,PlotLabel->StringForm["tstep=`1` t=`2`",t,t1[[t]]],PlotStyle->Blue],ListPlot[Flatten[Table[{ r1[[i]] ,plot2[[t,i]]} ,{i,1,rmax}],0],Joined->True, ImageSize->300,PlotRange->All,PlotStyle->Dark[Blue],PlotLabel->StringForm["tstep=`1` t=`2`",t,t1[[t]]]]},{t,plotmin,plotmax,plotint}],PlotLabel->Style[StringForm["fP=`1`, \[CapitalDelta]0=`2`,tmax=`3`,`4`-grid, `5`-range",fP,Delta,t1[[tmax]],grid,range],20]]*)

mylabel = 
  Style[Framed[
    StringForm[
     "fP=`1`, \[CapitalDelta]0=`2`, tmax=`3`, `4`-grid, `5`-range", 
     fP, Delta, Round[t1[[tmax]]], grid, range], 
    FrameStyle -> White, FrameMargins -> Tiny], 
   FontFamily -> "Arial", 21];
myep2 = Epilog -> {Text[
     Style["Re(\[Psi](r))", 18], 
     Scaled[{0.85, 0.85}]] };
myep1 = Epilog -> {Text[
     Style["|\[Psi](r)\!\(\*SuperscriptBox[\(|\), \(2\)]\)", 18], 
     Scaled[{0.80, 0.80}]] };

pairGridPlot[plot1_, plot2_] := 
 GraphicsGrid[
  Table[{ListPlot[
     Flatten[Table[{r1[[i]], plot1[[t, i]]}, {i, 1, rmax}], 0], 
     Joined -> True, ImageSize -> imgsize, PlotRange -> All, 
     PlotLabel -> 
      Style[StringForm["tstep=`1`, t=`2`", t, t1[[t]]], 23], 
     PlotStyle -> {Blue,Thick}, 
     AxesLabel -> {None,None}, 
     Ticks -> {{-10, -5, 0, 5, 10}, Automatic}, 
     TicksStyle -> {{FontSize -> 14, Black}, {FontSize -> 12, Black}},
      AxesOrigin -> {-12.5, 0}, GridLines -> {{0}, {0}}, 
     GridLinesStyle -> Directive[Black, Dashed],myep1], 
    ListPlot[
     Flatten[Table[{r1[[i]], plot2[[t, i]]}, {i, 1, rmax}], 0], 
     Joined -> True, ImageSize -> imgsize, PlotRange -> All, 
     PlotStyle -> {Dark[Teal],Thick}, 
     PlotLabel -> 
      Style[StringForm["tstep=`1`, t=`2`", t, t1[[t]]], 23], 
     AxesLabel -> {Style["r", FontSlant -> Italic, 18], 
       None},
      TicksStyle -> {{FontSize -> 12, Black}, {FontSize -> 12, 
        Black}}, AxesOrigin -> {-12.5, 0}, GridLines -> {{0}, {0}}, 
     GridLinesStyle -> Directive[Black, Dashed],myep2]}, {t, plotmin, 
    plotmax, plotint}], PlotLabel -> mylabel,Spacings -> {-30, 10}]


plot1=pairGridPlot[modatomsq1,reatom1];


outputname=StringTrim[filename,".h5"]<>".png" 

Export[outputname,plot1];
