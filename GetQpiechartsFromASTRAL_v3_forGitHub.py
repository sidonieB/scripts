############################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# This script plots the piecharts from an Astral tree
# I also created an easier alternative on R, no need for python or ete3 - see R scripts
#############################################################################################


### This script requires python and ete3 to be installed




### To prepare the input tree:

# generate Astral tree with quartet support. No branch lengths.

# if tree was generated with ASTRAL -t 2 option, need to do the following:

# replace ; by : (except the last ;)

# remove '

# add &&NHX: after each [

# This can be done in a text editor or in a unix terminal.

# invert annotations and branch lengths using the following command (use the command sed in a linux terminal - may not work exactly the same way in a mac - conservatively works with BL where unit is between 0 and 9, with 100 decimals and/or with 100 decimals and E-0 to E-999):

# sed 's/(\[[^\[]*\])(:[0-9]\.[0-9]{1,100}E{0,1}-{0,1}[0-9]{0,3})/\2\1/g' -r All_SpeciesTree_BP0_annotQ_rooted2.tre > All_SpeciesTree_BP0_annotQ_rooted3.tre



### To run the script on a remote server type "xvfb-run python GetQpiechartsFromASTRAL_v2.py InputFile.tre OutputFile.svg" in the terminal

### Locally "python GetQpiechartsFromASTRAL_v2.py InputFile.tre OutputFile.svg" should be enough

### You may first have to indicate the path to ete3, for instance by typing "export PATH=/home/user/anaconda_ete/bin:$PATH"

#############################################################################



### The script starts here:

from ete3 import Tree, TreeStyle, TextFace, AttrFace, NodeStyle, PieChartFace, faces
import sys


infile = sys.argv[1]
outfile = sys.argv[2]

t = Tree(infile)

###################################################################################################
#######ONLY THING TO CHANGE: OUTGROUP BUT YOU CAN ALSO ROOT THE TREE BEFORE########################
###################################################################################################

### root tree (change outgroup name as required, if rooting on a clade, use the commented get_common_ancestor function and replace the outgroup name by ancestor in the set_outgroup function)
#ancestor = t.get_common_ancestor("outtaxon1","outtaxon2","outtaxon3")

#t.set_outgroup("Corynocarpus_laevigata_AF206892-AF149001-AF479110-NC014807")

#############################################################
### Nothing more to change below this line EXCEPT if you are interested in making pies with posterior probabilities instead of quartet support, then see comments below.
#############################################################
# Set layout function
def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=10)
        faces.add_face_to_node(N, node, 0, position="aligned")

# Set tree style
ts = TreeStyle()
ts.show_leaf_name = False
ts.layout_fn = layout
nstyle = NodeStyle()
nstyle["shape"] = "sphere"
nstyle["size"] = 0
nstyle["fgcolor"] = "darkred"
ts.scale =  20


# Applies the same style to all nodes in the tree
for n in t.traverse():
   n.set_style(nstyle)

# Makes pie charts using the q values
#handle = open("output", "a")
for n in t.traverse():
   Q = []
   #P = []
   if hasattr(n,"q1"):
      Q.append(float(n.q1)*100)
   if hasattr(n,"q2"):
      Q.append(float(n.q2)*100)
   if hasattr(n,"q3"):
      Q.append(float(n.q3)*100)
      print n
      print Q
      Qtot = Q[0] + Q[1] + Q[2]
      print str(Qtot)
      qpie = PieChartFace(percents=Q, width=20, height=20, colors=["mediumblue", "darkorange", "dimgrey"], line_color=None)
      n.add_face(qpie, column=0, position = "branch-bottom")

### could do the same with p by commenting the section above and uncommenting the section below
   #if hasattr(n,"pp1"):
      #P.append(float(n.pp1)*100)
   #if hasattr(n,"pp2"):
   #   P.append(float(n.pp2)*100)
   #if hasattr(n,"pp3"):
   #   P.append(float(n.pp3)*100)
   #   handle.write(str(n) + "\n")
    #  handle.write(str(P[0]) + " ")
    #  handle.write(str(P[1]) + " ")
    #  handle.write(str(P[2]) + "\n")
    #  ppie = PieChartFace(percents=P, width=10, height=10, colors=["red", "cyan", "gray"], line_color=None)
    #  n.add_face(ppie, column=0, position = "branch-top")
     # n.dist = 1

#handle.close()

### ladderize the tree
t.ladderize(direction=1)

### export picture
t.render(outfile, w=183, units="mm", tree_style=ts)

#qpie gives the quartet support: percentage of quartets agreeing with the branch, with the second alternative RS|LO, and with the last alternative RO|LS.
#ppie gives local posterior probabilities, calculated using Q (see paper): one for the main topology OS|RL, and one for each of the two alternatives (RS|LO and RO|LS, in that order)

