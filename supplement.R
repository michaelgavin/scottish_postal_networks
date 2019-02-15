#### Supplemental script for "Infrastructural Semantics" ####
# Authors: Eric Gidal (Univeristy of Iowa), Michael Gavin (University of South Carolina)
# Corresponding author for all data and code: Gavin


           # # # # # # #     


#### SETUP AND INITIALIZATION ####

library(igraph)
library(empson)
# Set working directory
setwd("C:/Users/Michael Gavin/Desktop/scotland/")

##### Define functions #####
cos_sim = function(x, y, na.rm = T) {
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  (x %*% y) / (sqrt(x %*% x) * sqrt(y %*% y))
}

geodesicDistance <- function(long1, lat1, long2, lat2) {
  deg2rad <- function(deg) return(deg*pi/180)
  long1 = deg2rad(long1)
  lat1 = deg2rad(lat1)
  long2 = deg2rad(long2)
  lat2 = deg2rad(lat2)
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

map_words_export = function(words) {
  words = intersect(words, vocab)
  if (length(words) > 1) {
    vec = colSums(old[words,])
  } else {
    vec = old[words,]
  }
  freqs = vec / colSums(old)
  ids = names(freqs)
  old_ids = ids
  
  # New
  if (length(words) > 1) {
    vec = colSums(new[words,])
  } else {
    vec = new[words,]
  }
  freqs = c(freqs, vec / colSums(new))
  ids = names(vec)
  
  new_ids = ids
  
  df = paf[c(old_ids, new_ids),]
  
  df$ACCOUNT = "N"
  
  df$ACCOUNT[1:length(old_ids)] = "O"
  
  df$FREQ = freqs
  return(df)
}

footprint = function(words) {
  
  fp = list()
  fp$new = c()
  fp$old = c()
  
  lats = paf[place_ids,"Latitude"]
  lons = paf[place_ids,"Longitude"]
  
  if (length(words) > 1) {
    vec = colSums(new[words,])
  } else {
    vec = new[words,]
  }
  mean_lat = round(sum(vec * lats) / sum(vec), digits = 2)
  mean_lon = round(sum(vec * lons) / sum(vec), digits = 2)
  distances = c()
  for (i in 1:length(place_ids)) {
    d = geodesicDistance(mean_lon, mean_lat, lons[i], lats[i])
    distances = c(distances, d)
  }
  mean_d = round(mean(distances, na.rm = T), digits = 2)
  
  fp$new = c(sum(vec),mean_lat, mean_lon, mean_d)
  
  if (length(words) > 1) {
    vec = colSums(old[words,])
  } else {
    vec = old[words,]
  }
  mean_lat = round(sum(vec * lats) / sum(vec), digits = 2)
  mean_lon = round(sum(vec * lons) / sum(vec), digits = 2)
  distances = c()
  for (i in 1:length(place_ids)) {
    d = geodesicDistance(mean_lon, mean_lat, lons[i], lats[i])
    distances = c(distances, d)
  }
  mean_d = round(mean(distances, na.rm = T), digits = 2)
  
  fp$old = c(sum(vec),mean_lat, mean_lon, mean_d)
  
  return(fp)
  
}




##### LOAD DATA #####
# Load place authority file
paf = read.csv("paf.csv", stringsAsFactors = F)
paf = paf[!duplicated(paf$URL),]
paf = paf[!is.na(paf$Latitude),]
rownames(paf) = paf$URL

# Load edges
edges = read.csv("edges.csv", stringsAsFactors = F)
ord = unique(c(which(edges$Source %in% paf$URL == F),which(edges$Target %in% paf$URL == F)))
edges = edges[-ord,]

distances = c()
for (i in 1:nrow(edges)) {
  print(i)
  src = edges$Source[i]
  src_lat = paf$Latitude[which(paf$URL == src)][1]
  src_long = paf$Longitude[which(paf$URL == src)][1]
  
  targ = edges$Target[i]
  targ_lat = paf$Latitude[which(paf$URL == targ)][1]
  targ_long = paf$Longitude[which(paf$URL == targ)][1]
  
  d = geodesicDistance(src_long, src_lat, targ_long, targ_lat)
  d = round(d, digits = 2)
  distances = c(distances, d)
}

edges$Distance = distances


dmat = matrix(0, nrow(paf), nrow(paf))
rownames(dmat) = paf$URL
colnames(dmat) = paf$URL

for (i in 1:nrow(dmat)) {
  print(i)
  src_lat = paf$Latitude[i]
  src_long = paf$Longitude[i]
  for (j in 1:ncol(dmat)) {
    targ_lat = paf$Latitude[j]
    targ_long = paf$Longitude[j]
    d = geodesicDistance(src_long, src_lat, targ_long, targ_lat)
    d = round(d, digits = 2)
    dmat[i,j] = d
  }
}


# Load semantic models
load("new.rda")
load("old.rda")
new["hospital", ] = new["hospital",] + new["hofpital",]
old["hospital", ] = old["hospital",] + old["hofpital",]
new = new[-5967,]
old = old[-5967,]
stopwords = readLines("stopwords.txt")
place_ids = intersect(colnames(new), paf$URL)
vocab = setdiff(rownames(new), stopwords)
vocab = vocab[which(nchar(vocab) > 2)]

new = new[vocab, place_ids]
old = old[vocab, place_ids]

# Get basic stats for point distribution and perform regression analysis
e_dist = dmat["https://www.wikidata.org/wiki/Q23436", place_ids]
density = apply(dmat[place_ids, place_ids], 1, function(x) {length(which(x < 25)) / length(x)})
fit = lm(density ~ e_dist)



                   # # # # # # #     



#### Table 1 ####

selections = list()
selections[[1]] = c("Coach", "Omnibus")
selections[[2]] = "Carrier"
selections[[3]] = c("Mail", "Post")

ranges = list()
ranges[[1]] =  c(1794, 1799)
ranges[[2]] = c(1824,1820,1829,1830)
ranges[[3]] = c(1840)

network_data = matrix(0, 6, 9)
for (i in 1:3) {
  print(i)
  selection = selections[[i]]
  for (j in 1:3) {
    hits = which(edges$Date %in% ranges[[j]] & edges$Type %in% selection)
    df = edges[hits,]
    el = as.matrix(df[,1:2])
    hits = which(el[,2] != "" & el[,1] != "")
    el = el[hits,]
    g = graph.edgelist(el = el)
    column = j + ((i - 1) * 3)
    network_data[1,column] = vcount(g)
    network_data[2,column] = ecount(g)
    network_data[3,column] = sum(df$Distance, na.rm = T) / 1000
    network_data[4,column] = sum(as.integer(df$Travel_Frequecy), na.rm = T) / 1000
    network_data[5,column] = sum(df$Distance * as.integer(df$Travel_Frequecy), na.rm = T) / (1000 * 1000)
    if (vcount(g) > 0) {
      network_data[6,column] = centr_betw(g, directed = F)$centralization
    }
  }
}
network_data = round(network_data, digits = 2)
rownames(network_data) = c("location_count","route_count","total_distance","total_freq","capacity","centralization")
colnames(network_data) = c(paste("Coaches", c("1790s", "1820s", "1840s")),
                           paste("Carriers", c("1790s", "1820s", "1840s")),
                           paste("Mail", c("1790s", "1820s", "1840s")))



                      # # # # # # #
#### Figure 3 ####

# Data for Figure 3 can be loaded into R. The image itself was produced using QGIS
# from this data, then exported to SVG for editing/cleanup in Inkscape. Please
# note that the routes as listed in this data reflect the authors' judgments
# and were not directly noted in the eighteenth-century postal directories
figure3_data = read.csv("figure3_routes.csv", stringsAsFactors = F)

#### Figure 4 and Table 2 ####

selection = c(selections[[1]], selections[[3]])
par(mfrow = c(3,2))
for (i in 1:3) {
  hits = which(edges$Type %in% selection & edges$Date %in% ranges[[i]])
  df = edges[hits,]
  el = as.matrix(df[,1:2])
  hits = which(el[,2] != "" & el[,1] != "")
  el = el[hits,]
  g = graph.edgelist(el = el)
  
  geo_layout = as.matrix(paf[V(g)$name,c(5,4)])
  
  outsiders = which(geo_layout[,1] > 5 | geo_layout[,1] < -8 | geo_layout[,2] < 54)
  g = delete_vertices(g, V(g)[names(outsiders)])
  
  
  geo_layout = as.matrix(paf[V(g)$name,c(5,4)])
  
  V(g)$size = 2
  V(g)$label = ""
  V(g)$color = "darkred"
  
  plot(g, layout = geo_layout, edge.arrow.width = 0, vertex.frame.color = "darkred")
  plot(g, edge.arrow.width = 0, vertex.frame.color = "darkred")
}


# Table 2

load("degree_data.rda")
load("betweenness_data.rda")
d = as.matrix(degree_data[,3:30])
rownames(d) = paf[rownames(d),"Label"]
b = as.matrix(btw_data[,3:30])
rownames(b) = paf[rownames(b),"Label"]

# COLUMN KEY:
# 1. All dates. Full network.
# 2. All dates. "Coach" "Omnibus"
# 3. All dates. "Ship"
# 4. All dates. "Canal"
# 5. All dates. "Mail" "Post"
# 6. All dates. "Carrier"
# 7. All dates. "Rail"
# 8. 1790s. Full network.
# 9. 1790s. "Coach" "Omnibus"
# 10. 1790s. "Ship"
# 11. 1790s. "Canal"
# 12. 1790s. "Mail" "Post"
# 13. 1790s. "Carrier"
# 14. 1790s. "Rail"
# 15. 1820s. Full network.
# 16. 1820s. "Coach" "Omnibus"
# 17. 1820s. "Ship"
# 18. 1820s. "Canal"
# 19. 1820s. "Mail" "Post"
# 20. 1820s. "Carrier"
# 21. 1820s. "Rail"
# 22. 1840. Full network.
# 23. 1840. "Coach" "Omnibus"
# 24. 1840. "Ship"
# 25. 1840. "Canal"
# 26. 1840. "Mail" "Post"
# 27. 1840. "Carrier"
# 28. 1840. "Rail"

table2 = matrix("", 15, 6)
vec = sort(d[,1], decreasing = T)[1:15]
table2[,1] = names(vec)
table2[,2] = vec
vec = sort(d[,1] / (d[,8] + 1)^3, decreasing = T)[1:15]
table2[,3] = names(vec)
table2[,4] = as.integer(vec)
vec = sort(b[,1] / d[,1], decreasing = T)[1:15]
table2[,5] = names(vec)
table2[,6] = as.integer(vec)


                      # # # # # # #

#### Figures 5 through 8 ####

# For the semantic maps (figures 5-8) output was saved and exported to QGIS
# for visualization

# Create output to QGIS. These were each exported individually. In QGIS, point
# features were slightly scaled in size by term frequency to ease visualization
terms_list = list()
terms_list[[1]] = c("herring","cod","fishery","fisheries")
terms_list[[2]] = c("linen","cotton","flax","spinning","spindles","yarn","loom","looms")
terms_list[[3]] = c("mail","carrier","post","coach","coaches","carriers")
terms_list[[4]] = c("farmers","tenants","horses","potatoes","oats","pease","crops","harvest")

df = map_words_export(words = terms_list[[1]])
ord = which(df$ACCOUNT == "N" & df$FREQ > 0)
hits = which(df$FREQ[ord] > mean(df$FREQ[ord]))
write.csv(df[ord[hits],], file = "output_new.csv")

ord = which(df$ACCOUNT == "O" & df$FREQ > 0)
hits = which(df$FREQ[ord] > mean(df$FREQ[ord]))
write.csv(df[ord[hits],], file = "output_old.csv")

#### Table 3 -- Keywords by Category ####

####  Create List of Network Keywords  #### 

# Please note that the functions below perform several analyses not
# reported in the article.

# Get basic stats for point distribution and perform regression analysis
e_dist = dmat["https://www.wikidata.org/wiki/Q23436", place_ids]
density = apply(dmat[place_ids, place_ids], 1, function(x) {length(which(x < 25)) / length(x)})
fit = lm(density ~ e_dist)


# Define semantic similarity matrix
sim_mat = matrix(0, length(place_ids), length(place_ids))
rownames(sim_mat) = place_ids
colnames(sim_mat) = place_ids
for (i in 1:nrow(sim_mat)) {
  print(i)
  sim = empson::similarity(new, place_ids[i], margin = 2, fullResults = T)
  sim_mat[i,] = sim
}


# Semantic change of keywords
mat = new + old

# Calculate conceptual work for the New SAS
word_sims = c()
for (i in 1:nrow(new)) {
  word_sims = c(word_sims, cos_sim(new[i,], mat[i,]))
}
names(word_sims) = vocab
work = rowSums(new[,]) * (1 - word_sims)
sort(work, decreasing = T)[1:20]

# then for the old
word_sims = c()
for (i in 1:nrow(new)) {
  word_sims = c(word_sims, cos_sim(old[i,], mat[i,]))
}
names(word_sims) = vocab
work = rowSums(old) * (1 - word_sims)
sort(work, decreasing = T)[1:20]

# then check places
sims = c()
for (j in 1:ncol(new)) {
  sims = c(sims, cos_sim(new[,j], mat[,j]))
}
names(sims) = place_ids
work = colSums(new) * (1 - sims)

network_keywords = list()
# Get full network
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
d = degree(g)
d = d[names(d) %in% place_ids]
vec[names(d)] = d
network_keywords[[1]] = list()
names(network_keywords)[[1]] = "Full Network"
top_words = empson::similarity(new, vec)
new_sims = empson::similarity(new, vec, fullResults = T)
network_keywords[[1]]$new = top_words
top_words = empson::similarity(old, vec)
old_sims = empson::similarity(old, vec, fullResults = T)
network_keywords[[1]]$old = top_words


# Get local hubs
vec = rep(0, ncol(new))
names(vec) = colnames(new)
hits = which(edges$Date %in% c(1794, 1799))
el = edges[hits,]
g = graph.edgelist(el = as.matrix(el[,1:2]), directed = F)
d = degree(g)
d = d[names(d) %in% place_ids]
vec[names(d)] = d
vec1 = vec

vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
d = degree(g)
d = d[names(d) %in% place_ids]
vec[names(d)] = d
vec2 = vec

ord = which(vec1 < 3 & vec2 > 10)
nodes = names(ord)
vec = rep(0, ncol(new))
names(vec) = colnames(new)
vec[nodes] = vec2[nodes]
network_keywords[[2]] = list()
names(network_keywords)[[2]] = "Local Hubs"
top_words = empson::similarity(new, vec)
new_sims = empson::similarity(new, vec, fullResults = T)
network_keywords[[2]]$new = top_words
top_words = empson::similarity(old, vec)
old_sims = empson::similarity(old, vec, fullResults = T)
network_keywords[[2]]$old = top_words


# Get hubs as bridge nodes:
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
d = degree(g)
d = d[names(d) %in% place_ids]
btw = betweenness(g)
btw = btw[names(btw) %in% names(vec)]
vec[names(btw)] = btw / d
network_keywords[[3]] = list()
names(network_keywords)[[3]] = "Bridge Nodes"
top_words = empson::similarity(new, vec)
new_sims = empson::similarity(new, vec, fullResults = T)
network_keywords[[3]]$new = top_words
top_words = empson::similarity(old, vec)
old_sims = empson::similarity(old, vec, fullResults = T)
network_keywords[[3]]$old = top_words

# NOTE: For cities and regions below, the comparison vector not weighted by network stats -- all binary

# Cities
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
d = degree(g)
d = d[names(d) %in% place_ids]
d = d[d > (mean(d) * 2)]
vec[names(d)] = 1    # treated as binary vector
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[4]] = list()
names(network_keywords)[[4]] = "Cities"
network_keywords[[4]]$new = new_sims
network_keywords[[4]]$old = old_sims


# Five regions

# Central Belt
vec = fit$residuals
vec[vec < 0] = 0
vec[e_dist > mean(e_dist)] = 0
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[5]] = list()
names(network_keywords)[[5]] = "Central Belt"
network_keywords[[5]]$new = new_sims
network_keywords[[5]]$old = old_sims

# Aberdeenshire
vec = fit$residuals
vec[vec < 0] = 0
vec[which(e_dist > (max(e_dist) / 2) & e_dist < mean(e_dist))] = 0 # FOR ABERDEENSHIRE
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[6]] = list()
names(network_keywords)[[6]] = "Aberdeenshire"
network_keywords[[6]]$new = new_sims
network_keywords[[6]]$old = old_sims

# Outer Islands
vec = fit$residuals
vec[vec < 0] = 0
vec[e_dist < (max(e_dist) / 2)] = 0 
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[7]] = list()
names(network_keywords)[[7]] = "Outer Islands"
network_keywords[[7]]$new = new_sims
network_keywords[[7]]$old = old_sims

# Inner Countryside
vec = fit$residuals
vec[vec > 0] = 0
vec[e_dist > mean(e_dist)] = 0
vec = -1 * vec
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[8]] = list()
names(network_keywords)[[8]] = "Inner Countryside"
network_keywords[[8]]$new = new_sims
network_keywords[[8]]$old = old_sims

# Highlands
vec = fit$residuals
vec[vec > 0] = 0
vec[e_dist < mean(e_dist)] = 0
vec = -1 * vec
new_sims = empson::similarity(new, vec)
old_sims = empson::similarity(old, vec)
network_keywords[[9]] = list()
names(network_keywords)[[9]] = "Highlands"
network_keywords[[9]]$new = new_sims
network_keywords[[9]]$old = old_sims

# Convert to Keyword_DF
WORD = c()
VALUE = c()
CATEGORY = c()
SAS = c()
for (i in 1:length(network_keywords)) {
  x = network_keywords[[i]]
  for (j in 1:length(x)) {
    y = x[[j]]
    for (k in 1:length(y)) {
      WORD = c(WORD,names(y)[k])
      VALUE = c(VALUE, as.numeric(y[k]))
      CATEGORY = c(CATEGORY,names(network_keywords)[i])
      SAS = c(SAS, names(x)[j])
    }
  }
}
keywords_df = data.frame(CATEGORY, SAS, WORD, VALUE)
table3 = keywords_df[25:96,]

#### ADDITIONAL ANALYSES NOT REPORTED ####

# Measure semantic coherence of each region #
semantic_stats = list()

# Get local hubs
vec = rep(0, ncol(new))
names(vec) = colnames(new)
deg_dist = degree_list[[1]]
deg_dist = deg_dist[names(deg_dist) %in% names(vec)]
vec[names(deg_dist)] = deg_dist
vec1 = vec

vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
deg_dist = degree(g)
deg_dist = deg_dist[names(deg_dist) %in% names(vec)]
vec[names(deg_dist)] = deg_dist
vec2 = vec

ord = which(vec1 < 3 & vec2 > 10)
nodes = names(ord)
vec = rep(0, ncol(new))
names(vec) = colnames(new)
vec[nodes] = 1
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[1]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[1]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")


# Get hubs as bridge nodes:
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
deg_dist = degree(g)
deg_dist = deg_dist[names(deg_dist) %in% names(vec)]
btw = betweenness(g)
btw = btw[names(btw) %in% names(vec)]
vec[names(btw)] = btw / deg_dist
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[2]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[2]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Get full network
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
deg_dist = degree(g)
deg_dist = deg_dist[names(deg_dist) %in% names(vec)]
vec[names(deg_dist)] = deg_dist
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[3]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[3]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Central Belt
vec = fit$residuals
vec[vec < 0] = 0
vec[e_dist > mean(e_dist)] = 0
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[4]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[4]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Aberdeenshire
vec = fit$residuals
vec[vec < 0] = 0
vec[c(which(e_dist > (max(e_dist) / 2)), which(e_dist < mean(e_dist)))] = 0 # FOR ABERDEENSHIRE
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[5]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[5]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Outer Islands
vec = fit$residuals
vec[vec < 0] = 0
vec[e_dist < (max(e_dist) / 2)] = 0 
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[6]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[6]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Inner Countryside
vec = fit$residuals
vec[vec > 0] = 0
vec[e_dist > mean(e_dist)] = 0
vec = -1 * vec
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[7]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[7]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Highlands
vec = fit$residuals
vec[vec > 0] = 0
vec[e_dist < mean(e_dist)] = 0
vec = -1 * vec
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[8]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),                        
                        mean(1 - sims[ids]))
names(semantic_stats[[8]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")

# Cities
vec = rep(0, ncol(new))
names(vec) = colnames(new)
g = graph.edgelist(el = as.matrix(edges[,1:2]), directed = F)
deg_dist = degree(g)
deg_dist = deg_dist[names(deg_dist) %in% names(vec)]
deg_dist = deg_dist[deg_dist > (mean(deg_dist) * 2)]
vec[names(deg_dist)] = 1
vec = vec[vec > 0]
ids = names(vec)

self_sim = c()
other_sim = c()
for (i in 1:length(ids)) {
  print(i)
  id = ids[i]
  self_sim = c(self_sim,mean(sim_mat[id,ids]))
  other_sim = c(other_sim,mean(sim_mat[id, which(colnames(sim_mat) %in% ids == F)]))
}
vec = self_sim - other_sim
semantic_stats[[9]] = c(as.integer(length(vec)),
                        length(vec[vec > 0]) / length(vec),
                        mean(vec),
                        mean(self_sim),
                        mean(other_sim),
                        mean(work[ids]),
                        mean(1 - sims[ids]))
names(semantic_stats[[9]]) = c("N","Percent Self Sim", "Mean Difference", "Self Sim", "Other Sim","Work","Sem Distance")


names(semantic_stats) = c("Hubs","Bridges","Full", "Central Belt", "Aberdeenshire", "Outer Islands", "Inner Countryside", "Highlands", "Cities")

