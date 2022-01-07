myvarName = "Catch"
myvarDenomExpr = expression(slot(out_i,"RefPoint")$MSY[,,-(1:out_i@nyears),drop=FALSE])
ylabel = expression(mean~~Catch/MSY)
FUN = mean_noNA
scenNames=scenNames_user
OMNames= OMNames
MPNamesSub = MPNamesSub_user
MPNamesSub_legtext=MPNamesSub_user_legtext
legtext_o = 1 # order for plotting legend text
MSEoutNames = gsub("OM","MSE",OMNames)
colLabs = str_replace_all(gsub("OM_","",OMNames),"(?<=[a-z])(?=[A-Z])"," ")
rowLabs=scenNames_key[scenNames_user]

mar=c(1,1,1,1)
cols = setNames(c("red","black","deepskyblue2","deepskyblue2"),MPNamesSub)
ltys = setNames(c(1,1,1,2),MPNamesSub)
lwds = setNames(c(3,3,3,3),MPNamesSub)
types = setNames(c("l","l","l","l"),MPNamesSub)
xlabel = "Management year"
xlim = NULL
ylim=c(0,2)
assessmentYears = seq(10,40,by=10)
x.at = NULL
y.at = NULL
# leg.x=NULL
# leg.y =NULL
legend_panel=NULL
legend_cex= 1
legend_x="bottomleft"
legend_inset=c(0,-0.6)
label_axes="all" # "all" labels all panels, "leftbottom" labels y axes on the left and x axes on the bottom panels
hline=1 # values to pass into abline(h=hline)