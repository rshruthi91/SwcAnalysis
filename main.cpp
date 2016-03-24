#include <QString>
#include <QDebug>
#include <QFile>

#include <string>
#include <iostream>

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

  qDebug() << num_branches << "branches are there in this file";

  return EXIT_SUCCESS;
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
              unsigned int j = parentNodes.indexOf(curr_val);
              parentNodes.replace(j,-1);
              branch_count++;
            }
        }
    }
  return branch_count;
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
