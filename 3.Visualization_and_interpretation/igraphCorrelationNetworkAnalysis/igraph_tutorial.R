#这里是根据https://kateto.net/networks-r-igraph教程进行学习

# 创建网络图（create networks）
library(igraph)
g1 <- graph(edges = c(1,2, 
                      2,3, 
                      3,1), n = 3, directed = F)
plot(g1)
#上面这行代码生成了一个拥有三个边的无向图（undirected graph with three edges）函数中的数字是节点（vertex）的编号，边是1-2、2-3、3-1。

#检查一下g1在R中的存储类型
class(g1)

#在R中网络数据被专门存储为一个数据类型igraph，再看看g1在R中的存储方式
g1

#如果将节点数改成10个，方向（directed）为默认
g2 <- graph(edges = c(1,2, 
                      2,3, 
                      3,1), n = 10)
g2
plot(g2)

#也可以将节点的数字换成名字，而这时节点数也是不需要的
g3 <- graph(c("John","Jim",  
              "Jim","Jill",  
              "Jill","John"))
g3
plot(g3)

#在节点为名字的时候，我们可以只提供一些名字表示他们是独立的，并没有和其他点建立联系
g4 <- graph(c("John","Jim",  
              "Jim","Jack",  
              "Jim","Jack",
              "John","John"),
            isolates = c("Jesse","Janis","jennifer","Justin"))
plot(g4, 
     edge.arrow.size = .5,
     vertex.color = "gold",
     vertex.size = 15,
     vertex.frame.color = "grey",
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.label.dist = 2,
     edge.curved = 0.2
)

#小规模的网络图可以由一些描述（description）来生成：
#表述没有方向性的连接（undirected tie）
plot(graph_from_literal(a---b, b---c))

#+- 或者 -+ 表示指向左或右的方向性连接（directed tie）
plot(graph_from_literal(a--+b, b+--v))

#++ 表示双向的连接（symmetric tie）
plot(graph_from_literal(a+-+b, b+-+c))

#表示节点集合（sets of vertices）
plot(graph_from_literal(a:b:c---c:d:e))

plot(graph_from_literal(a-b-c-d-e-f, a-g-h-b, h-e:f:i, j))

#边、节点和网络属性（Edge, vertex, and network attributes）
g4 <- graph(c("John","Jim",  
              "Jim","Jack",  
              "Jim","Jack",
              "John","John"),
            isolates = c("Jesse","Janis","jennifer","Justin"))
#访问（assess）节点和边
E(g4)#访问边
V(g4)#访问节点

#直接查看网络矩阵（network matrix）
g4[]

#在网络图中的索引和矩阵类似
g4[1, ]#访问第一行

#对网络图、节点和边添加属性
V(g4)$name#name在我们创建网络图的时候以及自动生成了

V(g4)$gender <- c("male", "male", "male", 
                  "male", "female", "female", "male")

E(g4)$type <- "email"#边属性，对所有的边赋值（assign）“email”

E(g4)$weight <- 10#边权重，设置所有存在的边的权重为10

#查看属性
edge_attr(g4)#边属性

vertex_attr(g4)#节点属性

graph_attr(g4)#网络图属性

#其他的一些添加属性的方式
g4 <- set_graph_attr(g4, "name", "Email Network")
g4 <- set_graph_attr(g4, "something", "A thing")
graph_attr_names(g4)

graph_attr(g4)

#同样我们也可以删除这些属性
g4 <- delete_graph_attr(g4, "something")
graph_attr(g4)

#带着这些属性我们再绘制网络图
plot(g4,
     edge.arrow.size = .2, 
     vertex.label.color = "black",
     vertex.label.dist = 1.5,
     vertex.color = c("pink", "skyblue")[1+(V(g4)$gender=="male")],
)
#在上面这行代码中vertex.color=c(“pink”,“skyblue”)[1+(V(g4)($)gender==“male”)]：根据节点的性别（男/女）设置节点颜色为粉色/天蓝色。如果节点是男性，则颜色为天蓝色；否则颜色为粉色。这里使用了条件语句来选择颜色。其中，V(g4)(美元符号)gender表示g4图形对象中每个节点的gender属性，1表示该节点是男性，0表示该节点是女性。因为[1+(V(g4)(美元符号)gender==“male”)]返回的是一个向量，向量的值为1或2，所以用c(“pink”,“skyblue”)[1+(V(g4)(美元符号)gender==“male”)] 可以选择颜色。如果节点是男性，则返回 skyblue，否则返回 pink。

c("pink", "skyblue")[1 + (V(g4)$gender=="male")]

#从图中我们可以看到，在g4网络图中从Jim到Jack有两个边，还有一个环（loop）从John出发指向他自己。我们可以去掉从相同节点出发的环（loop）和多重边（multiple edges），同时我们也可以使用edge.attr.comb参数来指出边属性是怎样被一些选项给结合起来的，这些选项包括了sum,mean,prod (product 产出),min,max,first/last(选择第一或最后一个边的属性).“ignore”选项指出属性应该被忽略或丢弃
g4s <- simplify(g4,
                remove.multiple = T,
                remove.loops = F,
                edge.attr.comb = c(weight = "sum", type = "ignore"))

plot(g4s, vertex.label.dist = 1.5, edge.arrow.size = .2)
g4s

# 接下来让我们来看看上面结果对g4s这一网络图对象的描述：
# 关于igraph对象的描述是以四个字母开始的：
# D 或者 U，表明是有向或无向的网络图；
# N 表明这是一个节点有名字属性的命名图；
# W 表明这是一个边有权重属性的权重图；
# B 表明这是一个节点有类别属性的双边图（bipartite graph）。
# 下面的两个数字（7 3）分别表示的是在图中节点和边的数量。
# 然后还有关于节点和边的描述：
# (g/c) 表明图形水平（graph-level）是字符（character）属性
# (v/c) 表明节点水平（vertex-level）是字符（character）属性
# (e/n) 表明边水平（edge-level）是数字（numeric）属性

#简单星图（simple star graph）
st <- make_star(40)
plot(st, vertex.size = 10, vertex.label = NA)

#树图（tree graph）
tr <- make_tree(40, children = 3, mode = "undirected")
plot(tr, vertex.size = 10, vertex.label = NA)

#Watts-Strogatz 小世界模型 Watts-Strogatz小世界模型是一种用于生成小世界网络的模型。在该模型中，网络开始为一个具有N个节点的正则网格，每个节点与其k个最近邻节点相连。然后，以概率p（称为重连概率）重新连接节点的边，使得该节点与随机选择的另一个节点相连。这种重新连接的过程导致网络的局部结构（即节点周围的连接）保持不变，但是具有少量随机连接，从而增加了网络的全局连通性。当p很小时，网络仍然具有正则结构，而当p很大时，网络变成了随机网络。但是，在某个p值附近，网络具有小世界特性，即在保持局部结构的同时，具有短平均路径和高聚集性。
sw <- sample_smallworld(dim = 2, size = 10, nei = 1, p = 0.1)
plot(sw, 
     vertex.size = 6, 
     vertex.label = NA, 
     layout = layout_in_circle)

#Barabasi-Albert优先连接模型(用于生成无标度网络)(Barabasi-Albert preferential attachment model for scale-free graphs)Barabasi-Albert模型是一种用于生成无标度网络的模型，其核心思想是”优先连接”。在该模型中，网络开始只有几个节点，然后在每一步中，新节点以概率与现有节点连接，连接的概率与该节点的度数成正比。这意味着，具有更高度数的节点更有可能获得新的连接。这种优先连接的过程导致网络的度数分布遵循幂律分布，即存在少数高度连接的节点，大部分节点只连接少数其他节点。因此，该模型生成的网络通常被称为无标度网络，因为它们不遵循随机图模型中的泊松分布，而是具有一个长尾分布。这种模型在复杂网络的研究中具有广泛的应用。
ba <- sample_pa(n = 100, power = 1, m = 1, directed = F)
plot(ba, vertex.size = 6, vertex.label = NA)

#igraph can also give you some notable historical graphs.
zach <- graph("Zachary") # the Zachary carate club

plot(zach, vertex.size=10, vertex.label=NA)


#3. Reading network data from files
#3.1 DATASET 1: edgelist
nodes <- read.csv("netscix2016/Dataset1-Media-Example-NODES.csv", header=T, as.is=T)

links <- read.csv("netscix2016/Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)

head(nodes)

head(links)

nrow(nodes); length(unique(nodes$id))

nrow(links); nrow(unique(links[,c("from", "to")]))

#Notice that there are more links than unique from-to combinations. That means we have cases in the data where there are multiple links between the same two nodes. We will collapse all links of the same type between the same two nodes by summing their weights, using aggregate() by “from”, “to”, & “type”. We don’t use simplify() here so as not to collapse different link types.

links <- aggregate(links[,3], links[,-3], sum)

links <- links[order(links$from, links$to),]

colnames(links)[4] <- "weight"

rownames(links) <- NULL

#3.2 DATASET 2: matrix
#Two-mode or bipartite graphs have two different types of actors and links that go across, but not within each type. Our second media example is a network of that kind, examining links between news sources and their consumers.
nodes2 <- read.csv("netscix2016/Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)

links2 <- read.csv("netscix2016/Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)

head(nodes2)

head(links2)

#We can see that links2 is an adjacency matrix for a two-mode network:
links2 <- as.matrix(links2)

dim(links2)

dim(nodes2)


#4. Turning networks into igraph objects
#We start by converting the raw data to an igraph network object. Here we use igraph’s graph.data.frame function, which takes two data frames: d and vertices.

#d describes the edges of the network. Its first two columns are the IDs of the source and the target node for each edge. The following columns are edge attributes (weight, type, label, or anything else).
#vertices starts with a column of node IDs. Any following columns are interpreted as node attributes.
#4.1 Dataset 1
library(igraph)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 

class(net)

net

#We also have easy access to nodes, edges, and their attributes with:
E(net)       # The edges of the "net" object

V(net)       # The vertices of the "net" object

E(net)$type  # Edge attribute "type"

V(net)$media # Vertex attribute "media"

#Now that we have our igraph network object, let’s make a first attempt to plot it.
plot(net, edge.arrow.size=.4,vertex.label=NA)

#That doesn’t look very good. Let’s start fixing things by removing the loops in the graph.
net <- simplify(net, remove.multiple = F, remove.loops = T) 

#You might notice that we could have used simplify to combine multiple edges by summing their weights with a command like simplify(net, edge.attr.comb=list(weight="sum","ignore")). The problem is that this would also combine multiple edge types (in our data: “hyperlinks” and “mentions”).

#If you need them, you can extract an edge list or a matrix from igraph networks.
as_edgelist(net, names=T)

as_adjacency_matrix(net, attr="weight")

#Or data frames describing nodes and edges:
as_data_frame(net, what="edges")

as_data_frame(net, what="vertices")

#4.2 Dataset 2
#As we have seen above, this time the edges of the network are in a matrix format. We can read those into a graph object using graph_from_incidence_matrix(). In igraph, bipartite networks have a node attribute called type that is FALSE (or 0) for vertices in one mode and TRUE (or 1) for those in the other mode.
head(nodes2)
head(links2)
net2 <- graph_from_incidence_matrix(links2)

table(V(net2)$type)

#To transform a one-mode network matrix into an igraph object, use instead graph_from_adjacency_matrix().

#We can also easily generate bipartite projections for the two-mode network: (co-memberships are easy to calculate by multiplying the network matrix by its transposed matrix, or using igraph’s bipartite.projection() function).
net2.bp <- bipartite.projection(net2)

#We can calculate the projections manually as well:
as_incidence_matrix(net2)  %*% t(as_incidence_matrix(net2)) 

t(as_incidence_matrix(net2)) %*%   as_incidence_matrix(net2)

plot(net2.bp$proj1, vertex.label.color="black", vertex.label.dist=1,
     
     vertex.size=7, vertex.label=nodes2$media[!is.na(nodes2$media.type)])

plot(net2.bp$proj2, vertex.label.color="black", vertex.label.dist=1,
     
     vertex.size=7, vertex.label=nodes2$media[ is.na(nodes2$media.type)])


#5. Plotting networks with igraph
#5.1 Plotting parameters
# NODES	 
# vertex.color	 Node color
# vertex.frame.color	 Node border color
# vertex.shape	 One of “none”, “circle”, “square”, “csquare”, “rectangle”
#  “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
# vertex.size	 Size of the node (default is 15)
# vertex.size2	 The second size of the node (e.g. for a rectangle)
# vertex.label	 Character vector used to label the nodes
# vertex.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# vertex.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# vertex.label.cex	 Font size (multiplication factor, device-dependent)
# vertex.label.dist	 Distance between the label and the vertex
# vertex.label.degree	 The position of the label in relation to the vertex,
#  where 0 right, “pi” is left, “pi/2” is below, and “-pi/2” is above
# EDGES	 
# edge.color	 Edge color
# edge.width	 Edge width, defaults to 1
# edge.arrow.size	 Arrow size, defaults to 1
# edge.arrow.width	 Arrow width, defaults to 1
# edge.lty	 Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”,
#  3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
# edge.label	 Character vector used to label edges
# edge.label.family	 Font family of the label (e.g.“Times”, “Helvetica”)
# edge.label.font	 Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
# edge.label.cex	 Font size for edge labels
# edge.curved	 Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
# arrow.mode	 Vector specifying whether edges should have arrows,
#  possible values: 0 no arrow, 1 back, 2 forward, 3 both
# OTHER	 
# margin	 Empty space margins around the plot, vector with length 4
# frame	 if TRUE, the plot will be framed
# main	 If set, adds a title to the plot
# sub	 If set, adds a subtitle to the plo

#We can set the node & edge options in two ways - the first one is to specify them in the plot() function, as we are doing below.
# Plot with curved edges (edge.curved=.1) and reduce arrow size:

plot(net, edge.arrow.size=.4, edge.curved=.1)

# Set edge color to gray, and the node color to orange. 

# Replace the vertex label with the node names stored in "media"

plot(net, edge.arrow.size=.2, edge.curved=0,
     
     vertex.color="orange", vertex.frame.color="#555555",
     
     vertex.label=V(net)$media, vertex.label.color="black",
     
     vertex.label.cex=.7) 

#The second way to set attributes is to add them to the igraph object. Let’s say we want to color our network nodes based on type of media, and size them based on audience size (larger audience -> larger node). We will also change the width of the edges based on their weight.
# Generate colors based on media type:
colrs <- c("gray50", "tomato", "gold")
V(net)$color <- colrs[V(net)$media.type]
# Set node size based on audience size:
V(net)$size <- V(net)$audience.size*0.7
# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"
V(net)$label <- NA
# Set edge width based on weight:
E(net)$width <- E(net)$weight/6
#change arrow size and edge color:
E(net)$arrow.size <- .2
E(net)$edge.color <- "gray80"
E(net)$width <- 1+E(net)$weight/12

plot(net)

#We can also override the attributes explicitly in the plot:
plot(net, edge.color="orange", vertex.color="gray50") 

#It helps to add a legend explaining the meaning of the colors we used:
plot(net) 

legend(x=-1.5, y=-1.1, c("Newspaper","Television", "Online News"), pch=21,
       
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)

#Sometimes, especially with semantic networks, we may be interested in plotting only the labels of the nodes:
plot(net, vertex.shape="none", vertex.label=V(net)$media, 
     
     vertex.label.font=2, vertex.label.color="gray40",
     
     vertex.label.cex=.7, edge.color="gray85")

#Let’s color the edges of the graph based on their source node color. We can get the starting node for each edge with the ends() igraph function.
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- V(net)$color[edge.start]
plot(net, edge.color=edge.col, edge.curved=.1)  

#5.2 Network layouts
#Network layouts are simply algorithms that return coordinates for each node in a network.

#For the purposes of exploring layouts, we will generate a slightly larger 80-node graph. We use the sample_pa() function which generates a simple graph starting from one node and adding more nodes and links based on a preset level of preferential attachment (Barabasi-Albert model).
net.bg <- sample_pa(80) 
V(net.bg)$size <- 8
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- "" 
E(net.bg)$arrow.mode <- 0
plot(net.bg)

#You can set the layout in the plot function:
plot(net.bg, layout=layout_randomly)

#Or you can calculate the vertex coordinates in advance:
l <- layout_in_circle(net.bg)

plot(net.bg, layout=l)

#l is simply a matrix of x, y coordinates (N x 2) for the N nodes in the graph. You can easily generate your own:
l <- cbind(1:vcount(net.bg), c(1, vcount(net.bg):2))

plot(net.bg, layout=l)

#This layout is just an example and not very helpful - thankfully igraph has a number of built-in layouts, including:
# Randomly placed vertices
l <- layout_randomly(net.bg)

plot(net.bg, layout=l)

# Circle layout
l <- layout_in_circle(net.bg)

plot(net.bg, layout=l)

# 3D sphere layout
l <- layout_on_sphere(net.bg)

plot(net.bg, layout=l)

#Fruchterman-Reingold is one of the most used force-directed layout algorithms out there.

#Force-directed layouts try to get a nice-looking graph where edges are similar in length and cross each other as little as possible. They simulate the graph as a physical system. Nodes are electrically charged particles that repulse each other when they get too close. The edges act as springs that attract connected nodes closer together. As a result, nodes are evenly distributed through the chart area, and the layout is intuitive in that nodes which share more connections are closer to each other. The disadvantage of these algorithms is that they are rather slow and therefore less often used in graphs larger than ~1000 vertices. You can set the “weight” parameter which increases the attraction forces among nodes connected by heavier edges.
l <- layout_with_fr(net.bg)

plot(net.bg, layout=l)

#You will notice that the layout is not deterministic - different runs will result in slightly different configurations. Saving the layout in l allows us to get the exact same result multiple times, which can be helpful if you want to plot the time evolution of a graph, or different relationships – and want nodes to stay in the same place in multiple plots.
par(mfrow=c(2,2), mar=c(0,0,0,0))   # plot four figures - 2 rows, 2 columns

plot(net.bg, layout=layout_with_fr)

plot(net.bg, layout=layout_with_fr)

plot(net.bg, layout=l)

plot(net.bg, layout=l)

#By default, the coordinates of the plots are rescaled to the [-1,1] interval for both x and y. You can change that with the parameter rescale=FALSE and rescale your plot manually by multiplying the coordinates by a scalar. You can use norm_coords to normalize the plot with the boundaries you want.
l <- layout_with_fr(net.bg)

l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

par(mfrow=c(2,2), mar=c(0,0,0,0))

plot(net.bg, rescale=F, layout=l*0.4)

plot(net.bg, rescale=F, layout=l*0.6)

plot(net.bg, rescale=F, layout=l*0.8)

plot(net.bg, rescale=F, layout=l*1.0)

#Another popular force-directed algorithm that produces nice results for connected graphs is Kamada Kawai. Like Fruchterman Reingold, it attempts to minimize the energy in a spring system.
l <- layout_with_kk(net.bg)

plot(net.bg, layout=l)

#The LGL algorithm is meant for large, connected graphs. Here you can also specify a root: a node that will be placed in the middle of the layout.
plot(net.bg, layout=layout_with_lgl)

#Let’s take a look at all available layouts in igraph:
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# Remove layouts that do not apply to our graph.

layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]

par(mfrow=c(3,3), mar=c(1,1,1,1))

for (layout in layouts) {
  
  print(layout)
  
  l <- do.call(layout, list(net)) 
  
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }

#5.3 Improving network plots
#Notice that our network plot is still not too helpful. We can identify the type and size of nodes, but cannot see much about the structure since the links we’re examining are so dense. One way to approach this is to see if we can sparsify the network, keeping only the most important ties and discarding the rest.
hist(links$weight)

mean(links$weight)

sd(links$weight)

#There are more sophisticated ways to extract the key edges, but for the purposes of this exercise we’ll only keep ones that have weight higher than the mean for the network. In igraph, we can delete edges using delete_edges(net, edges):
cut.off <- mean(links$weight) 

net.sp <- delete_edges(net, E(net)[weight<cut.off])

plot(net.sp) 

#Another way to think about this is to plot the two tie types (hyperlink & mention) separately.

E(net)$width <- 1.5

plot(net, edge.color=c("dark red", "slategrey")[(E(net)$type=="hyperlink")+1],
     
     vertex.color="gray40", layout=layout.circle)

net.m <- net - E(net)[E(net)$type=="hyperlink"] # another way to delete edges

net.h <- net - E(net)[E(net)$type=="mention"]

# Plot the two links separately:

par(mfrow=c(1,2))

plot(net.h, vertex.color="orange", main="Tie: Hyperlink")

plot(net.m, vertex.color="lightsteelblue2", main="Tie: Mention")

# Make sure the nodes stay in place in both plots:
l <- layout_with_fr(net)

plot(net.h, vertex.color="orange", layout=l, main="Tie: Hyperlink")

plot(net.m, vertex.color="lightsteelblue2", layout=l, main="Tie: Mention")

#5.4 Interactive plotting with tkplot
#R and igraph allow for interactive plotting of networks. This might be a useful option for you if you want to tweak slightly the layout of a small graph. After adjusting the layout manually, you can get the coordinates of the nodes and use them for other plots.
tkid <- tkplot(net) #tkid is the id of the tkplot that will open

l <- tkplot.getcoords(tkid) # grab the coordinates from tkplot

tk_close(tkid, window.close = T)

plot(net, layout=l)

#5.5 Other ways to represent a network
#At this point it might be useful to provide a quick reminder that there are many ways to represent a network not limited to a hairball plot.
#For example, here is a quick heatmap of the network matrix:
netm <- get.adjacency(net, attr="weight", sparse=F)

colnames(netm) <- V(net)$media

rownames(netm) <- V(net)$media

palf <- colorRampPalette(c("gold", "dark orange")) 

heatmap(netm[,17:1], Rowv = NA, Colv = NA, col = palf(100), 
        
        scale="none", margins=c(10,10) )

#5.6 Plotting two-mode networks with igraph
#As with one-mode networks, we can modify the network object to include the visual properties that will be used by default when plotting the network. Notice that this time we will also change the shape of the nodes - media outlets will be squares, and their users will be circles.
V(net2)$color <- c("steel blue", "orange")[V(net2)$type+1]

V(net2)$shape <- c("square", "circle")[V(net2)$type+1]

V(net2)$label <- ""

V(net2)$label[V(net2)$type==F] <- nodes2$media[V(net2)$type==F] 

V(net2)$label.cex=.4

V(net2)$label.font=2

plot(net2, vertex.label.color="white", vertex.size=(2-V(net2)$type)*8) 

#Igraph also has a special layout for bipartite networks (though it doesn’t always work great, and you might be better off generating your own two-mode layout).
plot(net2, vertex.label=NA, vertex.size=7, layout=layout_as_bipartite) 

#Using text as nodes may be helpful at times:
plot(net2, vertex.shape="none", vertex.label=nodes2$media,
     
     vertex.label.color=V(net2)$color, vertex.label.font=2.5, 
     
     vertex.label.cex=.6, edge.color="gray70",  edge.width=2)

#6. Network and node descriptives
#6.1 Density
#The proportion of present edges from all possible edges in the network.
edge_density(net, loops=F)

ecount(net)/(vcount(net)*(vcount(net)-1)) #for a directed network

# 6.2 Reciprocity
# The proportion of reciprocated ties (for a directed network).
reciprocity(net)

dyad_census(net) # Mutual, asymmetric, and nyll node pairs

2*dyad_census(net)$mut/ecount(net) # Calculating reciprocity

# 6.3 Transitivity
# global - ratio of triangles (direction disregarded) to connected triples.
# local - ratio of triangles to connected triples each vertex is part of.

transitivity(net, type="global")  # net is treated as an undirected network

transitivity(as.undirected(net, mode="collapse")) # same as above

transitivity(net, type="local")

triad_census(net) # for directed networks 

# 6.4 Diameter
# A network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network. In igraph, diameter() returns the distance, while get_diameter() returns the nodes along the first found path of that distance.
# 
# Note that edge weights are used by default, unless set to NA.
diameter(net, directed=F, weights=NA)
diameter(net, directed=F)
diam <- get_diameter(net, directed=T)

diam

#Note that get_diameter() returns a vertex sequence. Note though that when asked to behaved as a vector, a vertex sequence will produce the numeric indexes of the nodes in it. The same applies for edge sequences.
class(diam)
as.vector(diam)

#Color nodes along the diameter:
vcol <- rep("gray40", vcount(net))

vcol[diam] <- "gold"

ecol <- rep("gray80", ecount(net))

ecol[E(net, path=diam)] <- "orange" 

# E(net, path=diam) finds edges along a path, here 'diam'

plot(net, vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0)

#6.5 Node degrees
#The function degree() has a mode of in for in-degree, out for out-degree, and all or total for total degree.
deg <- degree(net, mode="all")

plot(net, vertex.size=deg*3)

hist(deg, breaks=1:vcount(net)-1, main="Histogram of node degree")

#6.6 Degree distribution
deg.dist <- degree_distribution(net, cumulative=T, mode="all")

plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      
      xlab="Degree", ylab="Cumulative Frequency")

# 6.7 Centrality & centralization
# Centrality functions (vertex level) and centralization functions (graph level). The centralization functions return res - vertex centrality, centralization, and theoretical_max - maximum centralization score for a graph of that size. The centrality function can run on a subset of nodes (set with the vids parameter). This is helpful for large graphs where calculating all centralities may be a resource-intensive and time-consuming task.
# 
# Degree (number of ties)
degree(net, mode="in")

centr_degree(net, mode="in", normalized=T)

#Closeness (centrality based on distance to others in the graph)
#Inverse of the node’s average geodesic distance to others in the network.

closeness(net, mode="all", weights=NA) 

centr_clo(net, mode="all", normalized=T) 

#Eigenvector (centrality proportional to the sum of connection centralities)
#Values of the first eigenvector of the graph matrix.

eigen_centrality(net, directed=T, weights=NA)

centr_eigen(net, directed=T, normalized=T) 

#Betweenness (centrality based on a broker position connecting others)
#Number of geodesics that pass through the node or the edge.

betweenness(net, directed=T, weights=NA)

edge_betweenness(net, directed=T, weights=NA)

centr_betw(net, directed=T, normalized=T)

#6.8 Hubs and authorities
#The hubs and authorities algorithm developed by Jon Kleinberg was initially used to examine web pages. Hubs were expected to contain catalogs with a large number of outgoing links; while authorities would get many incoming links from hubs, presumably because of their high-quality relevant information.

hs <- hub_score(net, weights=NA)$vector

as <- authority_score(net, weights=NA)$vector

par(mfrow=c(1,2))

plot(net, vertex.size=hs*50, main="Hubs")

plot(net, vertex.size=as*30, main="Authorities")


# 7. Distances and paths
#Average path length: the mean of the shortest distance between each pair of nodes in the network (in both directions for directed graphs).
mean_distance(net, directed=F)
mean_distance(net, directed=T)
#We can also find the length of all shortest paths in the graph:
distances(net) # with edge weights

distances(net, weights=NA) # ignore weights

# We can extract the distances to a node or set of nodes we are interested in. Here we will get the distance of every media from the New York Times.
dist.from.NYT <- distances(net, v=V(net)[media=="NY Times"], to=V(net), weights=NA)

# Set colors to plot the distances:

oranges <- colorRampPalette(c("dark red", "gold"))

col <- oranges(max(dist.from.NYT)+1)

col <- col[dist.from.NYT+1]

plot(net, vertex.color=col, vertex.label=dist.from.NYT, edge.arrow.size=.6, 
     
     vertex.label.color="white")

#We can also find the shortest path between specific nodes. Say here between MSNBC and the New York Post:
news.path <- shortest_paths(net, 
                            
                            from = V(net)[media=="MSNBC"], 
                            
                            to  = V(net)[media=="New York Post"],
                            
                            output = "both") # both path nodes and edges

# Generate edge color variable to plot the path:

ecol <- rep("gray80", ecount(net))

ecol[unlist(news.path$epath)] <- "orange"

# Generate edge width variable to plot the path:

ew <- rep(2, ecount(net))

ew[unlist(news.path$epath)] <- 4

# Generate node color variable to plot the path:

vcol <- rep("gray40", vcount(net))

vcol[unlist(news.path$vpath)] <- "gold"

plot(net, vertex.color=vcol, edge.color=ecol, 
     
     edge.width=ew, edge.arrow.mode=0)

#Identify the edges going into or out of a vertex, for instance the WSJ. For a single node, use incident(), for multiple nodes use incident_edges()
inc.edges <- incident(net,  V(net)[media=="Wall Street Journal"], mode="all")

# Set colors to plot the selected edges.

ecol <- rep("gray80", ecount(net))

ecol[inc.edges] <- "orange"

vcol <- rep("grey40", vcount(net))

vcol[V(net)$media=="Wall Street Journal"] <- "gold"

plot(net, vertex.color=vcol, edge.color=ecol)


#We can also easily identify the immediate neighbors of a vertex, say WSJ. The neighbors function finds all nodes one step out from the focal actor.To find the neighbors for multiple nodes, use adjacent_vertices() instead of neighbors(). To find node neighborhoods going more than one step out, use function ego() with parameter order set to the number of steps out to go from the focal node(s).
neigh.nodes <- neighbors(net, V(net)[media=="Wall Street Journal"], mode="out")

# Set colors to plot the neighbors:

vcol[neigh.nodes] <- "#ff9d00"

plot(net, vertex.color=vcol)

# Special operators for the indexing of edge sequences: %–%, %->%, %<-%
# * E(network)[X %–% Y] selects edges between vertex sets X and Y, ignoring direction
# * E(network)[X %->% Y] selects edges from vertex sets X to vertex set Y
# * E(network)[X %->% Y] selects edges from vertex sets Y to vertex set X
# 
# For example, select edges from newspapers to online sources:

E(net)[ V(net)[type.label=="Newspaper"] %->% V(net)[type.label=="Online"] ]

#Co-citation (for a couple of nodes, how many shared nominations they have):
cocitation(net)

# 8. Subgroups and communities
# Before we start, we will make our network undirected. There are several ways to do that conversion:
# 
# We can create an undirected link between any pair of connected nodes (mode="collapse")
# Create undirected link for each directed one in the network, potentially ending up with a multiplex graph (mode="each")
# Create undirected link for each symmetric link in the graph (mode="mutual").
# In cases when we may have ties A -> B and B -> A ties collapsed into a single undirected link, we need to specify what to do with their edge attributes using the parameter ‘edge.attr.comb’ as we did earlier with simplify(). Here we have said that the ‘weight’ of the links should be summed, and all other edge attributes ignored and dropped.
net.sym <- as.undirected(net, mode= "collapse",
                         
                         edge.attr.comb=list(weight="sum", "ignore"))

#8.1 Cliques
#Find cliques (complete subgraphs of an undirected graph)
cliques(net.sym) # list of cliques       

sapply(cliques(net.sym), length) # clique sizes

largest_cliques(net.sym) # cliques with max number of nodes

vcol <- rep("grey80", vcount(net.sym))

vcol[unlist(largest_cliques(net.sym))] <- "gold"

plot(as.undirected(net.sym), vertex.label=V(net.sym)$name, vertex.color=vcol)

# 8.2 Community detection
# A number of algorithms aim to detect groups that consist of densely connected nodes with fewer connections across groups.
# 
# Community detection based on edge betweenness (Newman-Girvan)
# High-betweenness edges are removed sequentially (recalculating at each step) and the best partitioning of the network is selected. 
ceb <- cluster_edge_betweenness(net) 

dendPlot(ceb, mode="hclust")

plot(ceb, net) 

#Let’s examine the community detection igraph object:
class(ceb)
length(ceb)     # number of communities
membership(ceb) # community membership for each node
modularity(ceb) # how modular the graph partitioning is
crossing(ceb, net)   # boolean vector: TRUE for edges across communities

#High modularity for a partitioning reflects dense connections within communities and sparse connections across communities.

#Community detection based on based on propagating labels
# Assigns node labels, randomizes, than replaces each vertex’s label with the label that appears most frequently among neighbors. Those steps are repeated until each vertex has the most common label of its neighbors.
clp <- cluster_label_prop(net)

plot(clp, net)

#Community detection based on greedy optimization of modularity
cfg <- cluster_fast_greedy(as.undirected(net))

plot(cfg, as.undirected(net))

#We can also plot the communities without relying on their built-in plot:
V(net)$community <- cfg$membership

colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)

plot(net, vertex.color=colrs[V(net)$community])

#8.3 K-core decomposition
#The k-core is the maximal subgraph in which every node has degree of at least k. The result here gives the coreness of each vertex in the network. A node has coreness D if it belongs to a D-core but not to (D+1)-core.
kc <- coreness(net, mode="all")

plot(net, vertex.size=kc*6, vertex.label=kc, vertex.color=colrs[kc])

# 9. Assortativity and Homophily
# Homophily: the tendency of nodes to connect to others who are similar on some variable.
# 
# assortativity_nominal() is for categorical variables (labels)
# assortativity() is for ordinal and above variables
# assortativity_degree() checks assortativity in node degrees
assortativity_nominal(net, V(net)$media.type, directed=F)
assortativity(net, V(net)$audience.size, directed=F)
assortativity_degree(net, directed=F)
