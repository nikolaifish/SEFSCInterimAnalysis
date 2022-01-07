myvarName = "SB_SBMSY"
myvarDenomExpr = NULL
ylabel = myvarName
FUN = mean_noNA
OMNames=OMNames_user
scenNames=scenNames_user
# MSEoutNames = gsub("OM","MSE",OMNames)
MPNamesSub = c("SCA_10","SCA_1","iMP_avg_10","iMP_buffer_10")
MPNamesSub_legtext=MPNamesSub
legtext_o = c(2,3,4,1) # order for plotting legend text
MSEoutNames = paste0(rep(gsub("OM","MSE",OMNames),each=length(scenNames)),
                     "_",
                     rep(scenNames,length(OMNames))
)
colLabs = NULL
rowLabs = scenNames

mar=c(1,1,1,1)
cols = setNames(c("red","black","deepskyblue2","deepskyblue2"),MPNamesSub)
ltys = setNames(c(1,1,1,2),MPNamesSub)
lwds = setNames(c(3,3,3,3),MPNamesSub)
types = setNames(c("l","l","l","l"),MPNamesSub)
xlabel = "Management year"
xlim = NULL
ylim = NULL
assessmentYears = seq(10,40,by=10)
x.at = NULL
y.at = NULL
# leg.x=NULL
# leg.y =NULL
legend_panel=NULL
legend_cex= 1
legend_x="bottomleft"
legend_inset=c(0,-0.2)
label_axes="all" # "all" labels all panels, "leftbottom" labels y axes on the left and x axes on the bottom panels
hline=1 # values to pass into abline(h=hline)



ylabel = expression(mean~~SSB/SSBmsy)
FUN = mean
ylim=c(0,2)
MPNamesSub = MPNamesSub_user
MPNamesSub_legtext=MPNamesSub_user_legtext
rowLabs=scenNames_key[scenNames_user]