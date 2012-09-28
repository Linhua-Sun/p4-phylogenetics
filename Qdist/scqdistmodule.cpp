#include <Python.h>
#include "NewickParser.hpp"
#include "TreeUtil.hpp"
#include "QDist.hpp"
//#include <iostream>
 
int qdist(const char* input1, const char* input2) {


    Tree* tree1;
    Tree* tree2;

    NewickParser* parser = new NewickParser();
    tree1 = parser->Parse(input1);
    tree2 = parser->Parse(input2);

    TreeUtil::RenumberTreeAccordingToOther(tree2, tree1);    

    TreeUtil::CheckTree(tree1);
    TreeUtil::CheckTree(tree2);

    int result;
    result =  QDist::SubCubicQDist(tree1, tree2);
    //std::cout << "QDIST = " << result << std::endl;
    return result;

}
 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(scqdist)
{
    def("qdist", qdist);
}

