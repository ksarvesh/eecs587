#include <iostream>
#include <fstream>

using namespace std;



int main(int argc, char** argv){
    
    if(argc < 2){
        cout << "Usage: ./a.out <filename>" << endl;
        return(0);
    }

    string fileName(argv[1]);
    string fileName1= "../../test/" + fileName + ".txt";
    
    ifstream in(fileName.c_str());

    ofstream out(fileName1.c_str());
    
    string line;
    int numVertices, numEdges;
    
    if(!in.is_open()){
        cout<<"error"<<endl;
    }

    in>>numVertices;
    getline(in, line);
    out<<numVertices<<"\n";


    while(getline(in, line)){
        numEdges++;
    }

    out<<numEdges<<"\n";


    in.close();

    in.open(fileName.c_str());

    getline(in,line);


    while(getline(in,line)){
        out<<line<<"\n";
    }


    in.close();
    out.close();

    string fileName2= "../../base/" + fileName + ".txt";
    in.open(fileName.c_str());
    out.open(fileName2.c_str());

    getline(in, line);

    int from, to, capacity;

    while(getline(in, line)){
        sscanf(line.c_str(), "%d %d %d", &from, &to, &capacity);
        out<<"capacities["<<from<<"]["<<to<<"]="<<capacity<<";"<<"\n";
    }
    
    in.close();
    out.close();

    return(0);
}




