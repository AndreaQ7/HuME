library("ggspatial")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggrepel")
library("tmap")  

options(ggrepel.max.overlaps = Inf)
df=read.table("MappaItalia.txt",header=T)
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
    geom_sf(fill= "antiquewhite") +
    annotation_scale(location = "bl", width_hint = 0.5) +
    annotation_north_arrow(location = "bl", which_north = "true", 
        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
        style = north_arrow_fancy_orienteering) +
 geom_label_repel(data= df,aes(x=long, y=lat, label=Site),color = "grey22", fontface = "bold")+  geom_point(data=df,aes(x=long, y=lat))+ 
annotate(geom = "text", x = 18, y = 42, label = "Adriatic Sea", fontface = "italic", color = "#6699FF", size = 6) +
annotate(geom = "text", x = 12, y = 41, label = "Tyrrhenian Sea", fontface = "italic", color = "#6699FF", size = 6) +
coord_sf(xlim=c(9,20),ylim=c(35,44.5),expand=FALSE)+theme(panel.background = element_rect(fill = "aliceblue"))
