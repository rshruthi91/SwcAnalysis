#include <QString>
#include <QDebug>
#include <QFile>

#include <string>
#include <iostream>
#include <cmath>

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
int get_branches( QVector<int> parentNodes,QVector<int> *branch_nodes);
QVector<int> get_children(QVector<int> parent_nodes,int branch_val);
double absdiff(double n1, double n2);
void get_root_nodes(QVector<int> parentNodes,QVector<int> *root_nodes);
void get_terminal_nodes(QVector<vessel_node> vessels,QVector<int> parentNodes,QVector<int> *terminal_nodes);
void get_segment_seeds(QVector<int> parentNodes,QVector<int> *segment_seeds);
void calc_segment_area_vol(const vessel_node *node1, const vessel_node *node2, double *seg_vol, double *seg_lsa);


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
  int num_branches = get_branches(parent_nodes,&branch_nodes);
  if(num_branches < 0) return EXIT_FAILURE;

  qDebug() << num_branches << "branches are there in this trace";

  QVector<int> segment_seeds;
  segment_seeds = branch_nodes;
  get_segment_seeds(parent_nodes,&segment_seeds);
  int num_segments = segment_seeds.size() - 1;

  QVector<int> root_nodes;
  get_root_nodes(parent_nodes,&root_nodes);

  QVector<int> terminal_nodes;
  get_terminal_nodes(vessel_block.vessels,parent_nodes,&terminal_nodes);

  QVector<int> traversed_nodes;

//  QVector<int>::iterator vec_iter;
//  for(vec_iter = root_nodes.begin(); vec_iter != root_nodes.end();vec_iter++) {
//      if(terminal_nodes.indexOf(vec_iter) != -1) root_nodes.remove(vec_iter);
//    }

  double total_volume = 0.0;
  double total_surface_area = 0.0;
  for(int i=0; i < root_nodes.size();i++){
      bool end = false;
      int node1_num = root_nodes[i];
      int node2_num = 0; //Value '0' should never be used - dummy init.
      while(!end) {
          vessel_node node1 = vessel_block.vessels[node1_num];
          node2_num = parent_nodes.indexOf(node1_num) + 1;
          if( (branch_nodes.indexOf(node2_num) != -1) || (terminal_nodes.indexOf(node2_num) != -1) ){
              end = true;
              break;
            }
          vessel_node node2 = vessel_block.vessels[node2_num];

          double vol, lsa;
          calc_segment_area_vol(&node1,&node2,&vol,&lsa);

          total_volume+=vol;
          total_surface_area +=lsa;

          node1_num = node2_num;
        }
      traversed_nodes.append(root_nodes[i]);
    }

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

int get_branches( QVector<int> parentNodes,QVector<int> *branch_nodes){
  int curr_val,curr_cnt;
  int branch_count=0;
  for(int i=0;i<parentNodes.size();i++){
      curr_val =  parentNodes[i];
      if(curr_val != -1) {
          curr_cnt = parentNodes.count(curr_val);
          if(curr_cnt>4){
              qDebug() << "The SWC file is corrupt. More than 2 children for single parent";
              return -1;
            }
          if(curr_cnt>1) {
              branch_nodes->append(curr_val);
              parentNodes.replace(i,-1);
              //I can make a loop to replace one by one
              //but it should ot be necessary. >4 = error for now.
              unsigned int j = parentNodes.indexOf(curr_val);
              parentNodes.replace(j,-1);
              branch_count++;
            }
        }
    }
  return branch_count;
}

//Calculate the surface area and volume between two points
void calc_segment_area_vol(const vessel_node *node1, const vessel_node *node2, double *seg_vol, double *seg_lsa){
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
      QStringList somelist = line.split(' ');
      if(line[0] != '#'){
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
      else if(somelist.size() > 1){
          if(somelist[0] == "#XYZ"){
              vblock->x_dim = somelist[1].toInt();
              vblock->y_dim = somelist[2].toInt();
              vblock->z_dim = somelist[3].toInt();
            }
        }
    }
  file.close();
  return true;
}
