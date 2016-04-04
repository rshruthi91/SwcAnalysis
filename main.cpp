#include <QString>
#include <QDebug>
#include <QFile>
#include <QMap>

#include <string>
#include <iostream>
#include <cmath>
#include <limits>

# define PI 3.14159265358979323846  /* pi */

typedef struct {
   int  node_num;
   int type; // 0-Undefined, 1-Soma, 2-Axon, 3-Dendrite, 4-Apical_dendrite, 5-Fork_point, 6-End_point, 7-Custom
   double pos_x;
   double pos_y;
   double pos_z;
   double radius;
   int parent;
} vessel_node;

typedef struct{
 int x_dim;
 int y_dim;
 int z_dim;
 QVector<vessel_node>  vessels;
} VesselBlock;

bool read_swc_file(VesselBlock *vessel_blocks,QString Path);
int get_branches(QVector<int> parentNodes, QVector<int> *branch_nodes, QMap<int, QList<int> > *child_map);
QVector<int> get_children(QVector<int> parent_nodes,int branch_val);
double absdiff(double n1, double n2);
void get_root_nodes(QVector<int> parentNodes,QVector<int> *root_nodes);
void get_terminal_nodes(QVector<vessel_node> vessels,QVector<int> parentNodes,QVector<int> *terminal_nodes);
void get_segment_seeds(QVector<int> parentNodes,QVector<int> *segment_seeds);
void calc_segment_stats(const vessel_node *node1, const vessel_node *node2, double *seg_vol, double *seg_lsa, double *seg_len);
void copy_vnode(vessel_node in, vessel_node *out);

int main(int argc, char *argv[])
{
  VesselBlock vessel_block;

  std::string str;
  qDebug() << "Enter the name of swc file with path: " << endl;
  std::getline(std::cin, str);
  QString filename(str.c_str());
  if(!read_swc_file(&vessel_block,filename)) {
      std::cerr<< "Reading the SWC File Failed"<<std::endl;
      return EXIT_FAILURE;
    }
  qDebug() << "Reading the SWC File was successful";

  QVector<int> parent_nodes;
  for(int i=0; i < vessel_block.vessels.size();i++)
    parent_nodes.append(vessel_block.vessels[i].parent);

  QVector<int> branch_nodes;
  QMap<int, QList<int> > child_map;
  int num_branches = get_branches(parent_nodes,&branch_nodes, &child_map);
  if(num_branches < 0) return EXIT_FAILURE;

  QVector<int> segment_seeds;
  segment_seeds = branch_nodes;
  get_segment_seeds(parent_nodes,&segment_seeds);
  int num_segments = 0;

  QVector<int> root_nodes;
  get_root_nodes(parent_nodes,&root_nodes);

  QVector<int> terminal_nodes;
  get_terminal_nodes(vessel_block.vessels,parent_nodes,&terminal_nodes);
  segment_seeds+=terminal_nodes;

  QVector<int> traversed_nodes;

  double total_volume = 0.0;
  double total_surface_area = 0.0;
  double total_segment_length = 0.0;
  double max_segment_length = 0.0;
  double min_segment_length =  std::numeric_limits<double>::max();
  //FROM ROOTS TO ALL IMMEDIATE BRANCHES.
  //THIS PART IGNORES TERMINAL-ROOTS
  for(int i=0; i < root_nodes.size();i++){
      bool end = false;
      int node1_num = root_nodes[i];
      int node2_num = 0; //Value '0' should never be used - dummy init.
      vessel_node node1,node2;
      //If seed is a branch or terminal skip it
      if( (branch_nodes.indexOf(node1_num) >= 0) || (terminal_nodes.indexOf(node1_num) >= 0) ){
          end = true;
        }
      double vol=0.0;
      double lsa = 0.0;
      double seg_len = 0.0;
      while(!end) {
          num_segments++;
          copy_vnode(vessel_block.vessels[node1_num-1],&node1);
          node2_num = parent_nodes.indexOf(node1_num) + 1;
          if( (branch_nodes.indexOf(node1_num) >= 0) || (terminal_nodes.indexOf(node1_num) >= 0) ){
              end = true;
              continue;
            }
          copy_vnode(vessel_block.vessels[node2_num-1],&node2);

          calc_segment_stats(&node1,&node2,&vol,&lsa,&seg_len);

          total_volume+=vol;
          total_surface_area +=lsa;
          total_segment_length += seg_len;
          if(seg_len > max_segment_length) max_segment_length = seg_len;
          if(seg_len < min_segment_length) min_segment_length = seg_len;

          node1_num = node2_num;
        }
      traversed_nodes.append(root_nodes[i]);
    }
  qDebug() << endl;
  qDebug() << "Total Volume of segments from root to branch:" << total_volume << " voxel cubic units";
  qDebug() << "Total LSA of segments from root to branch:" << total_surface_area << "voxel square units";

  //BRANCH TO BRANCH TRACING
  foreach(int branch_node, branch_nodes){
      QList<int> children = child_map.value(branch_node);
      int node1_num=0,node2_num=0;
      vessel_node node1,node2;
      double vol=0.0;
      double lsa = 0.0;
      copy_vnode(vessel_block.vessels[branch_node-1],&node2);
      //Node1 is a branch, cannot be a terminal (else you wont know its a branch)
      foreach(int child_node, children){
          bool end = false;
          double seg_len = 0.0;
          node1_num = child_node;
          if( (branch_nodes.indexOf(node1_num) >= 0) || (terminal_nodes.indexOf(node1_num) >= 0) ){
              num_segments++;
              copy_vnode(vessel_block.vessels[child_node-1],&node1);
              calc_segment_stats(&node1,&node2,&vol,&lsa,&seg_len);
              total_volume+=vol;
              total_surface_area +=lsa;
              total_segment_length += seg_len;
              if(seg_len > max_segment_length) max_segment_length = seg_len;
              if(seg_len < min_segment_length) min_segment_length = seg_len;
              end = true;
              continue;
            }
          while(!end) {
              num_segments++;
              copy_vnode(vessel_block.vessels[node1_num-1],&node1);
              node2_num = parent_nodes.indexOf(node1_num) + 1;
              copy_vnode(vessel_block.vessels[node2_num-1],&node2);
              calc_segment_stats(&node1,&node2,&vol,&lsa,&seg_len);
              total_volume+=vol;
              total_surface_area +=lsa;
              total_segment_length += seg_len;
              if(seg_len > max_segment_length) max_segment_length = seg_len;
              if(seg_len < min_segment_length) min_segment_length = seg_len;
              if( (branch_nodes.indexOf(node2_num) >= 0) || (terminal_nodes.indexOf(node2_num) >= 0) ){
                  end = true;
                  continue;
                }
              node1_num = node2_num;
            }
        }
    }

  qDebug() << endl;
  qDebug() << "Total Volume of segments:" << total_volume << " voxel cubic units";
  qDebug() << "Total LSA of segments:" << total_surface_area << "voxel square units";

  double avg_seg_len = total_segment_length/num_segments;

  qDebug() << endl;
  qDebug() << "Avg Segment Length in structure" << avg_seg_len << "voxels";
  qDebug() << "Max Segment Length in structure" << max_segment_length << "voxels";
  qDebug() << "Min Segment Length in structure" << min_segment_length << "voxels";

  qDebug() << endl;
  qDebug() << "Total Number of Segments in Structure:" << num_segments;
  qDebug() << "Total Number of Branches in Structure:" << num_branches;
  qDebug() << "Total Number of Root Nodes in Structure:" << root_nodes.size();
  qDebug() << "Total Number of Terminals in Structure:" << terminal_nodes.size();

  return EXIT_SUCCESS;
}

QVector<int> get_children(QVector<int> parent_nodes,int branch_val){
  QVector<int> children;
  int from = 0;
  for(int i=0; i< parent_nodes.count(branch_val);i++){
      int inx = parent_nodes.indexOf(branch_val,from)+1;
      children.append(inx);
      from = inx+1;
    }
  return children;
}

void copy_vnode(vessel_node in, vessel_node *out){
  out->node_num = in.node_num;
  out->parent = in.parent;
  out->pos_x = in.pos_x;
  out->pos_y = in.pos_y;
  out->pos_z = in.pos_z;
  out->radius = in.radius;
  out->type = in.type;
}

double absdiff(double n1, double n2){
  if(n1>n2) return n1-n2;
  else return n2-n1;
}

void get_root_nodes(QVector<int> parentNodes,QVector<int> *root_nodes){
  for(int i =0; i<parentNodes.size();i++) {
      if(parentNodes[i] == -1) root_nodes->append(i+1);
    }
}

void get_terminal_nodes(QVector<vessel_node> vessels,QVector<int> parentNodes,QVector<int> *terminal_nodes){
  for(int i = 0; i < vessels.size(); i++){
      if(parentNodes.indexOf(vessels[i].node_num) == -1) terminal_nodes->append(vessels[i].node_num);
    }
}

void get_segment_seeds(QVector<int> parentNodes,QVector<int> *segment_seeds){
   for(int i=0;i<parentNodes.size();i++) { //Remember that i+1 is the node_num!!
       if(parentNodes[i] == -1) {
           segment_seeds->append(i+1);
         }
     }
}

int get_branches( QVector<int> parentNodes,QVector<int> *branch_nodes, QMap<int, QList<int> > *child_map){
  int curr_val,curr_cnt;
  int branch_count=0;
  for(int i=0;i<parentNodes.size();i++){
      curr_val =  parentNodes[i];
      if(curr_val != -1) {
          curr_cnt = parentNodes.count(curr_val);
          if(curr_cnt>1) {
              QList<int> child;
              branch_nodes->append(curr_val);
              child.append(i+1);
              parentNodes.replace(i,-1);
              if(curr_cnt > 2) {
                  for(int j = curr_cnt -2;j >0; j--) {
                      unsigned int jj = parentNodes.indexOf(curr_val);
                      parentNodes.replace(jj,-1);
                      child.append(jj+1);
                    }
                }
              child.append(parentNodes.indexOf(curr_val) + 1);
              child_map->insert(curr_val,child);
              child.clear();
              branch_count++;
            }
        }
    }
  return branch_count;
}

//Calculate the surface area and volume between two points
void calc_segment_stats(const vessel_node *node1, const vessel_node *node2, double *seg_vol, double *seg_lsa, double *seg_len){
  //Assume the vessel segment is a frustum of cone
  double x1_sq = (node1->pos_x)*(node1->pos_x);
  double y1_sq = (node1->pos_y)*(node1->pos_y);
  double z1_sq = (node1->pos_z)*(node1->pos_z);

  double x2_sq = (node2->pos_x)*(node2->pos_x);
  double y2_sq = (node2->pos_y)*(node2->pos_y);
  double z2_sq = (node2->pos_z)*(node2->pos_z);

  double h = std::sqrt(absdiff(x1_sq,x2_sq)+(y1_sq,y2_sq)+(z1_sq,z2_sq));
  double r1_sq = node1->radius*node1->radius;
  double r2_sq = node2->radius*node2->radius;
  double rdiff = absdiff(node1->radius,node2->radius);

  *seg_vol = PI*h*( r1_sq + r2_sq + node1->radius*node2->radius)/3;
  *seg_lsa = PI*(node1->radius+node2->radius)*sqrt(h*h + rdiff*rdiff);
  *seg_len = h;
}


bool read_swc_file(VesselBlock *vblock, QString inFileName) {

  QFile file(inFileName);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
      qDebug() << "Could not open the file in Read Only mode";
      return false;
    }
  QTextStream in(&file);

  while(!in.atEnd()){
      vessel_node node;
      QString line = in.readLine();
      QStringList somelist = line.split(' ', QString::SkipEmptyParts);
      if( (line[0] != '#') && (!line.isEmpty())){
          //qDebug() << somelist;
          node.node_num = somelist[0].toInt();
          node.type = somelist[1].toInt();
          node.pos_x = somelist[2].toDouble();
          node.pos_y = somelist[3].toDouble();
          node.pos_z = somelist[4].toDouble();
          node.radius = somelist[5].toDouble();
          node.parent = somelist[6].toInt();
          vblock->vessels.append(node);
         }
      /* IGNORE COMMENTS
      else if(somelist.size() > 1){
          if(somelist[0] == "#XYZ"){
              vblock->x_dim = somelist[1].toInt();
              vblock->y_dim = somelist[2].toInt();
              vblock->z_dim = somelist[3].toInt();
            }
        }*/
    }
  file.close();
  return true;
}
