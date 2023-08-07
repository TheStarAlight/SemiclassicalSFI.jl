using Colors
using ColorSchemes

r = cat(range(1.0,0.0,32),zeros(64),range(0.0,1.0,64),ones(64),range(1.0,0.5,32); dims=1)
g = cat(range(1.0,0.0,32),range(0.0,1.0,64),ones(64),range(1.0,0.0,64),zeros(32); dims=1)
b = cat(ones(96),range(1.0,0.0,64),zeros(96); dims=1)

RGB_list = map((r,g,b)->RGB(r,g,b), r,g,b)

loadcolorscheme(:ijet, RGB_list)