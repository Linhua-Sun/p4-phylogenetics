#include <iostream>
#include <pyublas/numpy.hpp>
#include <vector>
#include <boost/dynamic_bitset.hpp>

class BtKB
{
public:
  int nAllTax;
  boost::dynamic_bitset<> bK;
  int isInOther;

  BtKB(int nAllTax): nAllTax(nAllTax), isInOther(0) {
    boost::dynamic_bitset<> b(nAllTax);
    bK = b;
  }
};

class No
{
public:
  int nodeNum;
  int nAllTax;
  int taxNum;
  int isLeaf;
  boost::dynamic_bitset<> bitKey;
  //boost::dynamic_bitset<> bitKeyB;
  BtKB* bKK;
  No* parent;
  No* leftChild;
  No* sibling;
  int isRelevant;
  int isInOther;

  No(int nNum, int nAllTax): nodeNum(nNum), nAllTax(nAllTax), 
			     taxNum(-1), isLeaf(0), 
			     //bitKey(NULL), bKK(NULL),
			     parent(NULL), leftChild(NULL), sibling(NULL),
			     isRelevant(0), isInOther(0)
  {
    boost::dynamic_bitset<> bK(nAllTax);
    bitKey = bK;
    //bKK = new BtKB(nAllTax);
  }

  void setInternalBits() {
    bitKey = this->leftChild->bitKey;
    //std::cout << "  a No->setInternalBits() Setting node " << this->nodeNum << " to " << bitKey << std::endl;
    No* sib = this->leftChild->sibling;
    while(sib != NULL) {
      //std::cout << "  b sib node " << sib->nodeNum << ", bitKey=" << sib->bitKey << std::endl;
      bitKey |= sib->bitKey;
      //bitKey->set();
      sib = sib->sibling;
    }
    //std::cout << "  No->setInternalBits() Setting node " << this->nodeNum;
    //std::cout << " (isLeaf=" << this->isLeaf << ") to " << bitKey << std::endl;
  }

  void dump() {
      std::cout << this->nodeNum << " ";

      No* nn = this->parent;
      std::cout << "parent=";
      if(nn != NULL) {
  	std::cout << nn->nodeNum;
      } else {
  	std::cout << "None";
      }
      std::cout << "   ";

      nn = this->leftChild;
      std::cout << "leftChild=";
      if(nn != NULL) {
  	std::cout << nn->nodeNum;
      } else {
  	std::cout << "None";
      }
      std::cout << "   ";

      nn = this->sibling;
      std::cout << "sibling=";
      if(nn != NULL) {
  	std::cout << nn->nodeNum;
      } else {
  	std::cout << "None";
      }
      std::cout << "   ";
      
      std::cout << "isLeaf=";
      std::cout << this->isLeaf;
      std::cout << "   ";
    
      std::cout << "taxNum=";
      std::cout << this->taxNum;
      std::cout << "   ";
    
      std::cout << "bitKey=";
      std::cout << this->bitKey;
      std::cout << "   ";
    
      
      std::cout << std::endl;
  }
};

class Tr
{
public:
  int nNo;
  int nTax;
  pyublas::numpy_vector<int> postOrder;
  int nAllTax;
  double beta;
  //std::vector<boost::dynamic_bitset<>*> bitKeys; 
  std::vector<BtKB*> bitKeyBs; 
  boost::dynamic_bitset<> taxBits;
  std::size_t firstOnePos;
  std::vector<No*> nodes;

  Tr(int nN, int nTax, pyublas::numpy_vector<int> pO, int nAllTax, double beta): 
    nNo(nN), nTax(nTax), postOrder(pO), nAllTax(nAllTax), beta(beta), //taxBits(NULL), 
    firstOnePos(-1)
  {
    for(int i=0; i < nNo; i++) {
      No* newNo = new No(i, nAllTax);
      nodes.push_back(newNo);
    }
    //std::cout << "Made " << nodes.size() << " nodes." << std::endl;
    postOrder = pO;
    boost::dynamic_bitset<> tB(nAllTax);
    taxBits = tB;
  } 

  void setParent(int nNum, int oNum) {
    No* theNode=nodes[nNum];
    No* theOther=nodes[oNum];
    theNode->parent = theOther;
  }
  void setLeftChild(int nNum, int oNum) {
    No* theNode=nodes[nNum];
    No* theOther=nodes[oNum];
    theNode->leftChild = theOther;
  }
  void setSibling(int nNum, int oNum) {
    No* theNode=nodes[nNum];
    No* theOther=nodes[oNum];
    theNode->sibling = theOther;
  }

  void setNodeTaxNum(int nNum, int taxNum) {
    No* theNode=nodes[nNum];
    theNode->isLeaf = 1;
    theNode->taxNum = taxNum;
    theNode->bitKey[taxNum] = 1;
  }

  void setInternalBits() {
    No* n;
    for(int i = 0; i < nNo; i++) {
      n = nodes[this->postOrder[i]];
      if(n->isLeaf != 1) {
  	n->setInternalBits();
      }
    }
  }

  void setInTreeTaxBits() {
    No* n;
    for(int i = 0; i < nNo; i++) {
      n = nodes[i];
      if(n->isLeaf == 1) {
  	taxBits[n->taxNum] = 1;
      }
    }
    const std::size_t first_one = taxBits.find_first();
    firstOnePos = first_one;
  }

  void maybeFlipInTreeBits() {
    No* n;
    std::size_t fop;
    for(int i = 0; i < (nNo - 1); i++) {
      n = nodes[this->postOrder[i]];
      if(n->isLeaf == 0) {
  	fop = n->bitKey.find_first();
  	if(fop == this->firstOnePos) {
  	  //std::cout << "maybeFlipInTreeBits()  flipping node " << n->nodeNum << " " << n->bitKey << " fop=" << fop;
  	  n->bitKey.flip();
	  n->bitKey &= this->taxBits;
  	  //std::cout << " new bit key: " << n->bitKey << std::endl;
  	}
      }
    }
  }

  void wipePointers() {
    No* n;
    for(int i = 0; i < nNo; i++) {
      n = nodes[i];
      n->parent = NULL;
      n->leftChild = NULL;
      n->sibling = NULL;
    }
  }

  int flagRelevantSplits(Tr* iT) {
    // for bigT only
    No* n;
    std::size_t nOnes;
    int nRel;
    nRel = 0;
    for(int i = 0; i < (nNo - 1); i++) {
      n = nodes[postOrder[i]];
      if(n->isLeaf == 0) {
  	n->isRelevant = 0;
  	n->bKK->bK = n->bitKey;
  	if(n->bKK->bK[iT->firstOnePos] == 1) {
  	  n->bKK->bK.flip();
  	  //std::cout << "    flipping ..." << std::endl;
  	}
  	n->bKK->bK &= iT->taxBits;
  	//std::cout << "BigT node " << n->nodeNum << " bitKey " << n->bitKey << " iT->taxBits " << iT->taxBits;
  	//std::cout << " set bKK->bK " << n->bKK->bK << std::endl;
	
  	nOnes = n->bKK->bK.count();
  	if((nOnes >= 2ul) && (nOnes <= (unsigned long)(iT->nTax - 2))) {
  	  n->isRelevant = 1;
	  //std::cout << "     isRelevant" << std::endl;
  	  n->bKK->isInOther = 0;
  	  nRel += 1;
  	}
      }
    }
    return nRel;
  }

  void fillBitKeyBsVector() {
    No* n;
    int alreadyIn;
    BtKB* vBK;

    this->bitKeyBs.clear();
    for(int i = 0; i < (nNo - 1); i++) {
      n = nodes[postOrder[i]];
      if((n->isLeaf == 0) && (n->isRelevant == 1)) {
	alreadyIn = 0;
	for(std::size_t i2=0; i2 < this->bitKeyBs.size(); i2++) {
	  vBK = this->bitKeyBs[i2];
	  if(vBK->bK == n->bKK->bK) {
	    alreadyIn = 1;
	    break;
	  }
	}
	if(alreadyIn == 0) {
	  this->bitKeyBs.push_back(n->bKK);
	}
      }
    }
    //std::cout << "Tr.fillBitKeyBsVector(), got " << this->bitKeyBs.size() << " BtKB objects.\n";
  }

  double testIsInOther(Tr* iT) {
    No* itn;
    std::size_t i;
    BtKB* bKK;

    for(i = 0; i < this->bitKeyBs.size(); i++) {
      bKK = this->bitKeyBs[i];
      for(int i2 = 0; i2 < (iT->nNo - 1); i2++) {
	itn = iT->nodes[iT->postOrder[i2]];
	if((itn->isLeaf == 0) && (itn->isInOther == 0) && (bKK->bK == itn->bitKey)) {
	  //std::cout << "Bit key " << itn->bitKey << " is found in both" << std::endl;
	  bKK->isInOther = 1;
	  itn->isInOther = 1;
	  break;
	}
      }
      
    }

    double total;
    total = 0.0;
    for(i = 0; i < this->bitKeyBs.size(); i++) {
      bKK = this->bitKeyBs[i];
      if(bKK->isInOther == 0) {
  	total += 1.0;
      }
    }
    for(int i2 = 0; i2 < (iT->nNo - 1); i2++) {
      itn = iT->nodes[iT->postOrder[i2]];
      if((itn->isLeaf == 0) && (itn->isInOther == 0)) {
  	total += 1.0;
      }
    }

    return total * iT->beta;
  }

  void dump() {
    std::cout << "Tr dump().  nNo=" << this->nNo << ", firstOnePos="<< this->firstOnePos;
    std::cout << ", taxBits=" << this->taxBits << std::endl;
    No* n;
    for(int i=0; i< nNo; i++) {
      n = nodes[i];
      n->dump();
    }
  }


};

class Stm
{
public:
  int nTax;
  Tr* bigT;
  std::vector<Tr*> inTrees;
  std::vector<boost::dynamic_bitset<>*> singleBits;

  Stm(int nTax): nTax(nTax), 
		 bigT(NULL) 
  {
  } 
  
  Tr setBigT(int nN, int nTax, pyublas::numpy_vector<int> pO) {
    Tr* t;
    No* n;
    t = new Tr(nN, nTax, pO, this->nTax, 0.0);  // last arg, beta, does not matter for bigT.
    for(std::size_t i=0; i < (unsigned long)(t->nNo); i++) {
      n = t->nodes[i];
      n->bKK = new BtKB(nTax);
    }
    bigT = t;
    return *bigT;
  }

  Tr appendInTree(int nN, int nTax, pyublas::numpy_vector<int> pO, double beta) {
    Tr* inTree;
    inTree = new Tr(nN, nTax, pO, this->nTax, beta);
    this->inTrees.push_back(inTree);
    return *inTree;
  }

  void setInTreeInternalBits() {
    Tr* t;
    for(std::size_t i = 0; i < inTrees.size(); i++) {
      t = inTrees[i];
      //std::cout << "Stm->setInTreeInternalBits()  setting internal bits for inTree " << i << std::endl;
      t->setInternalBits();
    }
  }

  void setInTreeTaxBits() {
    Tr* t;
    for(std::size_t i = 0; i < inTrees.size(); i++) {
      t = inTrees[i];
      t->setInTreeTaxBits();
      //std::cout << "Stm.setInTreeTaxBits() setting taxBits for inTree " << i << " to " << t->taxBits;
      //std::cout << " firstOnePos " << t->firstOnePos << std::endl;
    }
  }

  void maybeFlipInTreeBits() {
    Tr* t;
    for(std::size_t i = 0; i < inTrees.size(); i++) {
      t = inTrees[i];
      t->maybeFlipInTreeBits();
    }
  }

  void setBigTInternalBits() {
    //std::cout << "Stm->setBigTInternalBits()  setting internal bits for bigT " << std::endl;
    this->bigT->setInternalBits();
  }

  void wipeBigTPointers() {
    //std::cout << "Stm->wipeBigTPointers()" << std::endl;
    this->bigT->wipePointers();
  }

  double getSymmDiff() {
    Tr* iT;
    int ret;
    double retd;
    No* n;
    int nNum;
    double sd;
    sd = 0.0;
    for(std::size_t i = 0; i < inTrees.size(); i++) {
      iT = inTrees[i];
      for(nNum=0; nNum < iT->nNo; nNum++) {
	n = iT->nodes[nNum];
	n->isInOther = 0;
      }
      for(nNum=0; nNum < this->bigT->nNo; nNum++) {
	n = this->bigT->nodes[nNum];
	n->bKK->isInOther = 0;
	n->isRelevant = 0;
      }
      ret = this->bigT->flagRelevantSplits(iT);
      //std::cout << "got ret=" << ret << " relevant" << std::endl;

      this->bigT->fillBitKeyBsVector();
      retd = this->bigT->testIsInOther(iT);
      //std::cout << "got this sd=" << retd << std::endl;
      sd += retd;
    }
    return sd;
  }

  void dump() {
    std::cout << "Stm dump(), nTax=" << this->nTax << std::endl;
    if(this->bigT == NULL) {
      std::cout << "  bigT is not set" << std::endl;
    } else {
      std::cout << "---- bigT dump ---------" << std::endl;
      this->bigT->dump();
      std::cout << "---- end of bigT dump --" << std::endl;
    }
    Tr* t;
    for(std::size_t i = 0; i < inTrees.size(); i++) {
      t = this->inTrees[i];
      std::cout << "---- inTree " << i << " dump ---------" << std::endl;
      t->dump();
      std::cout << "--------" << std::endl;
    }
  }

};


#include <boost/python.hpp> 
//#include <boost/python/module.hpp>
//#include <boost/python/def.hpp>
using namespace boost::python;
 
BOOST_PYTHON_MODULE(p4stm)
{

  class_<Tr>("Tr", init<int, int, pyublas::numpy_vector<int>, int, double >())
    .def("setParent", &Tr::setParent)
    .def("setLeftChild", &Tr::setLeftChild)
    .def("setSibling", &Tr::setSibling)
    .def("setNodeTaxNum", &Tr::setNodeTaxNum)
    .def("dump", &Tr::dump)
    ;
  class_<Stm>("Stm", init<int>())
    .def("setBigT", &Stm::setBigT)
    .def("appendInTree", &Stm::appendInTree)
    .def("setInTreeInternalBits", &Stm::setInTreeInternalBits)
    .def("setInTreeTaxBits", &Stm::setInTreeTaxBits)
    .def("maybeFlipInTreeBits", &Stm::maybeFlipInTreeBits)
    .def("setBigTInternalBits", &Stm::setBigTInternalBits)
    .def("wipeBigTPointers", &Stm::wipeBigTPointers)
    .def("getSymmDiff", &Stm::getSymmDiff)
    .def("dump", &Stm::dump)
    ;
}

