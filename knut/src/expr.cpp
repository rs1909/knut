#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <list>
#include <cctype>
#include <cmath>

#define KNUTSYS_H 

#include "knerror.h"
#include "matrix.h"
#include "expr.h"
#include "exprsystem.h"

namespace ExpTree
{

struct TokenPrecedence
{
  int precedence;
  bool leftAssoc;
  size_t nargs;
};

class Token
{
  public:
    Token (TokenType t, const std::string& n) : type(t), name(n) { }
    bool isId() const { return type == TokenType::Symbol; }
    bool isOp() const { return (type != TokenType::Invalid)&&(type != TokenType::LeftParen)&&(type != TokenType::RightParen)&&(type != TokenType::Symbol)&&(type != TokenType::Number)&&(type != TokenType::Comma); }
    bool isSym() const { return type == TokenType::Symbol; }
    bool isLeft() const { return type == TokenType::LeftParen; }
    bool isRight() const { return type == TokenType::RightParen; }
// variables
    TokenType type;
    std::string name;
};

// NodeNumber
bool NodeNumber::operator != (const Node& node) const 
{ 
  const NodeNumber* nptr = dynamic_cast<const NodeNumber*>(&node);
  if (nptr) return (value != nptr->value);
  else return true;
}

bool NodeNumber::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
    return value < static_cast<const NodeNumber* >(&node)->value;
  else return false;
}

// NodeSymbol
bool NodeSymbol::operator != (const Node& node) const 
{ 
  const NodeSymbol* nptr = dynamic_cast<const NodeSymbol*>(&node);
  if (nptr) return (name != nptr->name);
  else return true;
}

bool NodeSymbol::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
    return name < static_cast<const NodeSymbol* >(&node)->name;
  else return false;
}

void NodeSymbol::replaceSymbol(const Node& sym, const Node& node, Node** parent)
{
  if (*parent != this) P_MESSAGE1("Wrong parent!\n");
  const NodeSymbol* symref = dynamic_cast<const NodeSymbol*>(&sym);
  if (symref)
  {
    if (name == symref->name)
    {
      Node* cp = node.copy();
      (*parent)->deleteTree ();
      delete *parent;
      *parent = cp;
    }
  }
}

void NodeSymbol::replaceSymbol (const Node& node, Node** parent)
{
  if (node.type == TokenType::Equals)
  {
    const Node* left = node.getChild(0);
    const Node* right = node.getChild(1);
    if (left->type == TokenType::Symbol)
    {
      replaceSymbol (*left, *right, parent);
    }
  } else P_MESSAGE1("Not a definition\n");
}

// NodeVar

bool NodeVar::operator != (const Node& node) const
{ 
  const NodeVar* nptr = dynamic_cast<const NodeVar*>(&node);
  if (nptr) return (idx != nptr->idx);
  else return true;
}

bool NodeVar::operator < (const Node& node) const
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
  {
    if (idx < static_cast<const NodeVar* >(&node)->idx) return true;
    else return false;
  }
  else return false;
}

void NodeVar::replaceSymbol (const Node& sym, const Node& node, Node** parent)
{
  if (*parent != this) P_MESSAGE1("Wrong parent!\n");
  const NodeVar* symref = dynamic_cast<const NodeVar*>(&sym);
  if (symref)
  {
    if (idx == symref->idx)
    {
      Node* cp = node.copy ();
      (*parent)->deleteTree ();
      delete *parent;
      *parent = cp;
    }
  }
}

// NodePar

bool NodePar::operator != (const Node& node) const
{ 
  const NodePar* nptr = dynamic_cast<const NodePar*>(&node);
  if (nptr) return (idx != nptr->idx);
  else return true;
}
bool NodePar::operator < (const Node& node) const
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
    return idx < static_cast<const NodePar* >(&node)->idx;
  else return false;
}

void NodePar::replaceSymbol (const Node& sym, const Node& node, Node** parent)
{
  if (*parent != this) P_MESSAGE1("Wrong parent!\n");
  const NodePar* symref = dynamic_cast<const NodePar*>(&sym);
  if (symref)
  {
    if (idx == symref->idx)
    {
      Node* cp = node.copy ();
      (*parent)->deleteTree ();
      delete *parent;
      *parent = cp;
    }
  }
}

//NodeFunction
Node* NodeFunction::copy () const
{
  NodeFunction* cp = new NodeFunction(name, args.size());
  for(size_t k=0; k<args.size(); k++) cp->args[k] = args[k]->copy();
  return cp;
}

void NodeFunction::print (std::ostream& out) const 
{
  out << name << "(";
  for(size_t k=0; k<args.size(); k++) 
  {
    P_ERROR_X1( args[k] != nullptr, "Incerdibly wrong!");
    args[k]->print (out); 
    if (k<args.size()-1) out << ","; 
  }
  out << ")";
}

bool NodeFunction::operator != (const Node& node) const 
{ 
  const NodeFunction* nptr = dynamic_cast<const NodeFunction*>(&node);
  if (nptr)
  {
    if (args.size() != nptr->args.size()) return true;
    for (size_t k=0; k<args.size(); k++) if (*args[k] != *(nptr->args[k])) return true;
  }
  else return true;
  return false;
}

bool NodeFunction::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
  {
    if (name < static_cast<const NodeFunction* >(&node)->name) return true;
    if (name == static_cast<const NodeFunction* >(&node)->name)
    {
      for (size_t k=0; k<args.size(); k++)
      {
        if (*args[k] < *(static_cast<const NodeFunction* >(&node)->args[k])) return true;
        else if (*args[k] == *(static_cast<const NodeFunction* >(&node)->args[k])) continue;
        else return false;
      }
    } else return false;
  } 
  return false;
}

void NodeFunction::replaceSymbol(const Node& sym, const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
}

void NodeFunction::replaceSymbol(const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (node, &args[k]);
  if (node.type == TokenType::Equals)
  {
    const Node* left = node.getChild(0);
    const Node* right = node.getChild(1);
    if (left->type == TokenType::Function)
    {
      const NodeFunction* macro = static_cast<const NodeFunction* >(left);
      if ( (name == macro->name)&&(children() == macro->children()) )
      {
        Node* cp = right->copy ();
        for (size_t k=0; k<args.size(); k++)
        {
          if (macro->args[k]->type == TokenType::Symbol)
          {
            cp->replaceSymbol (*macro->args[k], *args[k], &cp);
          }
          else P_MESSAGE1("Not a symbol");
        }
        (*parent)->deleteTree ();
        delete *parent;
        *parent = cp;
      }
    }
  }
}

void NodeFunction::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
  if (name == nm) lst.push_back (this);
}

size_t NodeFunction::stackCount(size_t sp, size_t& max_stack) 
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> stackCount(sp+k, max_stack); }
  if (max_stack < spi) max_stack = spi; 
  return sp;
}

static double localSign (double v)
{
  if (v > 0.0) return 1.0;
  else if(v == 0.0) return 0.0;
  else return -1.0;
}

static double localHeaviside (double v)
{
  if (v > 0.0) return 1.0;
  else return 0.0;
}

// NodeFunctionA1
class NodeFunctionA1 : public NodeFunction
{
  public:
    NodeFunctionA1 (const std::string& nm);
    Node* copy () const;
    size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
    size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  private:
    std::function<double(double)> function;
};

NodeFunctionA1::NodeFunctionA1 (const std::string& nm) : NodeFunction (nm, 1)
{
  if (nm == "abs") function = fabs;
  else if (nm == "exp") function = exp;
  else if (nm == "log") function = log;
  else if (nm == "sin") function = sin;
  else if (nm == "cos") function = cos;
  else if (nm == "tan") function = tan;
  else if (nm == "asin") function = asin;
  else if (nm == "acos") function = acos;
  else if (nm == "atan") function = atan;
  else if (nm == "sinh") function = sinh;
  else if (nm == "cosh") function = cosh;
  else if (nm == "tanh") function = tanh;
  else if (nm == "asinh") function = asinh;
  else if (nm == "acosh") function = acosh;
  else if (nm == "atanh") function = atanh;
  else if (nm == "sign") function = localSign;
  else if (nm == "heaviside") function = localHeaviside;
}

Node* NodeFunctionA1::copy () const
{
  NodeFunctionA1* cp = new NodeFunctionA1(NodeFunction::name);
  for(size_t k=0; k<this->args.size(); k++) cp->args[k] = this->args[k]->copy();
  cp->function = function;
  return cp;
}

// NodeFunctionA2
class NodeFunctionA2 : public NodeFunction
{
  public:
    NodeFunctionA2 (const std::string& nm);
    Node* copy () const;
    size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const ;
    size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const;
  private:
    std::function<double(double,double)> function;
};

NodeFunctionA2::NodeFunctionA2 (const std::string& nm) : NodeFunction (nm, 2)
{
  if (nm == "pow") function = pow;
}

Node* NodeFunctionA2::copy () const
{
  NodeFunctionA2* cp = new NodeFunctionA2(NodeFunction::name);
  for(size_t k=0; k<this->args.size(); k++) cp->args[k] = this->args[k]->copy();
  cp->function = function;
  return cp;
}

// NodeAdd

Node* NodeAdd::copy () const
{
  NodeAdd* cp = new NodeAdd(args.size());
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->copy();
    cp->mul[k] = mul[k];
  }
  return cp;
}

void NodeAdd::print (std::ostream& out) const 
{
  P_ERROR_X1 (args.size() > 0, "terribly wrong");
  out << "(";
  for (size_t k=0; k<args.size(); k++)
  {
    P_ERROR_X1( args[k] != nullptr, "Incerdibly wrong!");
    if (mul[k] < 0.0) out << "-";
    else if (k != 0) out << "+";
    if ((mul[k] != 1.0) && (mul[k] != -1.0)) out << fabs(mul[k]) << "*";
    args[k]->print (out);
  }
  out << ")";
}

bool NodeAdd::operator != (const Node& node) const 
{ 
  const NodeAdd* nptr = dynamic_cast<const NodeAdd*>(&node);
  if (nptr)
  {
    if (args.size() != nptr->args.size()) return true;
    // the arguments can be in any order
    // see if one matches
    std::list<size_t> found;
    for (size_t k=0; k<args.size(); k++)
    {
      bool differs = true;
      for (size_t l=0; l<args.size(); l++) 
        if ( !((*args[k] != *(nptr->args[l])) || (mul[l] != nptr->mul[l])) ) differs = false;
      if (differs) return true;
    }
  }
  else return true;
  return false;
}

bool NodeAdd::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
  {
    if (args.size() < static_cast<const NodeAdd* >(&node)->args.size()) return true;
    else if (args.size() == static_cast<const NodeAdd* >(&node)->args.size())
    {
      for (size_t k=0; k<args.size(); k++)
      {
        if (mul[k] < static_cast<const NodeAdd* >(&node)->mul[k]) return true;
        else if (mul[k] == static_cast<const NodeAdd* >(&node)->mul[k]) 
        {  
          if (*args[k] < *(static_cast<const NodeAdd* >(&node)->args[k])) return true;
          else if (*args[k] == *(static_cast<const NodeAdd* >(&node)->args[k])) continue;
          else return false;
        } else return false;
      }
    } else return false;
  }
  return false;
}

void NodeAdd::replaceSymbol(const Node& sym, const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
}

void NodeAdd::replaceSymbol(const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++)
  {
    P_ERROR_X1 (args[k] != nullptr, "terribly wrong");
    args[k]->replaceSymbol (node, &args[k]);
  }
}

void NodeAdd::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
}

size_t NodeAdd::stackCount(size_t sp, size_t& max_stack) 
{
  size_t spi = sp + 1;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> stackCount(sp+k+1, max_stack); }
  if (max_stack < spi) max_stack = spi;
  return sp;
}

void NodeAdd::sort ()
{
  std::vector<size_t> order(args.size());
  for (size_t k=0; k<args.size(); k++) order[k] = k;
  std::sort (order.begin(), order.end(), [&](size_t a, size_t b) {
      return *(args[a]) < *(args[b]);
  });
  std::vector<Node*> tmp_a(args.size());
  std::vector<double> tmp_m(args.size());
  for (size_t k=0; k<args.size(); k++)
  {
    tmp_a[k] = args[order[k]];
    tmp_m[k] = mul[order[k]];
  }
  args = std::move(tmp_a);
  mul = std::move(tmp_m);
}

// NodeTimes

Node* NodeTimes::copy () const
{
  NodeTimes* cp = new NodeTimes(args.size());
  cp->smul = smul;
  for(size_t k=0; k<args.size(); k++)
  {
    P_ERROR_X1 (args[k] != nullptr, "terribly wrong");
    cp->args[k] = args[k] -> copy ();
    cp->divide[k] = divide[k];
  }
  return cp;
}

void NodeTimes::print (std::ostream& out) const 
{
  P_ERROR_X1 (args.size() > 0, "terribly wrong");
  out << "(";
  if (smul != 1.0)
  {
    if (smul < 0.0) out << "-";
    if (smul != -1.0) out << fabs(smul) << "*";
  }
  for(size_t k=0; k<args.size(); k++)
  {
    P_ERROR_X1 (args[k] != nullptr, "terribly wrong");
    if (k != 0) { if (divide[k]) out << "/"; else out << "*"; }
    else if (divide[k]) out << "1/";
    args[k]->print (out);
  }
  out << ")";
}

bool NodeTimes::operator != (const Node& node) const 
{ 
  const NodeTimes* nptr = dynamic_cast<const NodeTimes*>(&node);
  if (nptr)
  {
    if (args.size() != nptr->args.size()) return true;
    for(size_t k=0; k<args.size(); k++) if ((*args[k] != *(nptr->args[k])) || (divide[k] != nptr->divide[k])) return true;
  }
  else return true;
  return false;
}

bool NodeTimes::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
  {
    if (args.size() < static_cast<const NodeTimes* >(&node)->args.size()) return true;
    else if (args.size() == static_cast<const NodeTimes* >(&node)->args.size())
    {
      if (smul < static_cast<const NodeTimes* >(&node)->smul) return true;
      else if (smul == static_cast<const NodeTimes* >(&node)->smul)
      {
        for (size_t k=0; k<args.size(); k++)
        {
          if (divide[k] < static_cast<const NodeTimes* >(&node)->divide[k]) return true;
          else if (divide[k] == static_cast<const NodeTimes* >(&node)->divide[k]) 
          {  
            if (*args[k] < *(static_cast<const NodeTimes* >(&node)->args[k])) return true;
            else if (*args[k] == *(static_cast<const NodeTimes* >(&node)->args[k])) continue;
            else return false;
          } else return false;
        }
      } else return false;
    } else return false;
  }
  return false;
}

void NodeTimes::replaceSymbol(const Node& sym, const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
}

void NodeTimes::replaceSymbol(const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++)
  {
    P_ERROR_X4 (args[k] != nullptr, "terribly wrong: ", k, " of ", args.size() );
    args[k]->replaceSymbol (node, &args[k]);
  }
}

void NodeTimes::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
}

size_t NodeTimes::stackCount(size_t sp, size_t& max_stack) 
{
  size_t spi = sp + 1;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> stackCount(sp+k+1, max_stack); }
  if (max_stack < spi) max_stack = spi; 
  return sp;
}

void NodeTimes::sort ()
{
  std::vector<size_t> order(args.size());
  for (size_t k=0; k<args.size(); k++) order[k] = k;
  std::sort (order.begin(), order.end(), [&](size_t a, size_t b) {
      return *(args[a]) < *(args[b]);
  });
  std::vector<Node*> tmp_a(args.size());
  std::vector<bool> tmp_d(args.size());
  for (size_t k=0; k<args.size(); k++)
  {
    tmp_a[k] = args[order[k]];
    tmp_d[k] = divide[order[k]];
  }
  args = std::move(tmp_a);
  divide = std::move(tmp_d);
}

void NodeTimes::addArgument (size_t pos, Node* arg, double m, bool div)
{
  args[pos] = arg; 
  if (div) smul /= m; else smul *= m;
  divide[pos] = div;
}

// NodePower

Node* NodePower::copy () const
{
  NodePower* cp = new NodePower();
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->copy();
  }
  return cp;
}

void NodePower::print (std::ostream& out) const 
{
  out << "pow(";
  P_ERROR_X1( args[0] != nullptr, "Incerdibly wrong!");
  P_ERROR_X1( args[1] != nullptr, "Incerdibly wrong!");
  args[0]->print (out); out << ",";
  args[1]->print (out);
  out << ")";
}

bool NodePower::operator != (const Node& node) const 
{ 
  const NodePower* nptr = dynamic_cast<const NodePower*>(&node);
  if (nptr)
  {
    if (args.size() != nptr->args.size()) return true;
    for(size_t k=0; k<args.size(); k++) if (*args[k] != *(nptr->args[k])) return true;
  }
  else return true;
  return false;
}

bool NodePower::operator < (const Node& node) const 
{
  if (this->type < node.type) return true; 
  else if (this->type == node.type)
  {
    for (size_t k=0; k<args.size(); k++)
    {
      if (*args[k] < *(static_cast<const NodePower* >(&node)->args[k])) return true;
      else if (*args[k] == *(static_cast<const NodePower* >(&node)->args[k])) continue;
      else return false;
    }
  } 
  return false;
}

void NodePower::replaceSymbol(const Node& sym, const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
}

void NodePower::replaceSymbol(const Node& node, Node** parent)
{
  for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (node, &args[k]);
}

void NodePower::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
}

Node* NodePower::derivative (const Node* var) const 
{
  // D[a^b] = a^(b-1) * ( (D[a] * b) + (a * D[b] * log[a]))
  //          a  b -1      D[a]   b     a   D[b]      a
  //          |  \ /          \  /       \   |       /
  //          \   s1           t1         \  |    fun 
  //           \   /              \        \-t2--/
  //             pw                \         /
  //               \                \-- s2--/    
  //                \-----\             / 
  //                       \---  t3 ---/
  // root is t3
  const Node* a = args[0];
  const Node* b = args[1];
  NodeTimes* t1 = new NodeTimes(2);
  NodeTimes* t2 = new NodeTimes(3);
  NodeTimes* t3 = new NodeTimes(2);
  NodeAdd* s1 = new NodeAdd(2);
  NodeAdd* s2 = new NodeAdd(2);
  NodeFunctionA1* fun = new NodeFunctionA1("log");
  NodePower* pw = new NodePower;
  
  s1 -> addArgument (0, b -> copy (), 1.0);
  s1 -> addArgument (1, new NodeNumber(1.0), -1.0);
  pw -> args[0] = a -> copy();
  pw -> args[1] = s1;
  t1 -> addArgument (0, a -> derivative (var), 1.0, false);
  t1 -> addArgument (1, b -> copy (), 1.0, false);
  fun -> addArgument(0, a -> copy ());
  t2 -> addArgument (0, a -> copy (), 1.0, false);
  t2 -> addArgument (1, b -> derivative (var), 1.0, false);
  t2 -> addArgument (2, fun, 1.0, false);
  s2 -> addArgument (0, t1, 1.0);
  s2 -> addArgument (1, t2, 1.0);
  t3 -> addArgument (0, pw, 1.0, false);
  t3 -> addArgument (1, s2, 1.0, false);
  return t3;
}

void NodePower::optimize (Node** parent) 
{ 
  for (size_t k=0; k<args.size(); k++) args[k]->optimize (&args[k]);
  // if both of them are numbers ...
  NodeNumber* a0 = dynamic_cast<NodeNumber*>(args[0]);
  NodeNumber* a1 = dynamic_cast<NodeNumber*>(args[1]);
  if ((a0 != nullptr)&&(a1 != nullptr))
  {
    NodeNumber* res = new NodeNumber( pow (a0->getValue(), a1->getValue()) );
    (*parent)->deleteTree ();
    delete (*parent);
    (*parent) = res;
    return;
  }
  if (a1 != nullptr)
  {
    if (a1 -> getValue () == 1.0)
    {
      Node* tmp = args[0];
      delete args[1];
      delete (*parent);
      (*parent) = tmp;
      return;
    }
    if (a1 -> getValue () == 0.0)
    {
      NodeNumber* res = new NodeNumber (1.0);
      (*parent)->deleteTree ();
      delete (*parent);
      (*parent) = res;
      return;
    }
  }
}

size_t NodePower::stackCount(size_t sp, size_t& max_stack) 
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> stackCount(sp+k, max_stack); }
  if (max_stack < spi) max_stack = spi; 
  return sp;
}

class NodePlus : public Node
{
public:
  NodePlus () : Node(TokenType::UnaryPlus), arg(nullptr) {}
  void deleteTree () { arg->deleteTree (); delete arg; }
  Node* copy () const
  {
    NodePlus* cp = new NodePlus();
    cp->arg = arg->copy();
    return cp;
  }
  size_t children () const { return 1; }
  const Node* getChild (size_t n) const { return arg; }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const 
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const 
  {
    P_ERROR_X1( arg != nullptr, "Incerdibly wrong!");
    out << "(";
    out << "+";
    arg->print (out);
    out << ")";
  }
  bool operator != (const Node& node) const 
  { 
    const NodePlus* nptr = dynamic_cast<const NodePlus*>(&node);
    if (nptr) return (*arg != *(nptr->arg));
    else return true;
  }
  bool operator < (const Node& node) const 
  {
    if (this->type < node.type) return true; 
    else if (this->type == node.type)
    {
      if (*arg < *(static_cast<const NodePlus* >(&node)->arg)) return true;
      else return false;
    } 
    return false;
  }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent)
  {
    arg->replaceSymbol (sym, node, &arg);
  }
  void replaceSymbol (const Node& node, Node** parent)
  {
    arg->replaceSymbol (node, &arg);
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
  {
    arg->findFunction (lst, nm);
  }
  Node* derivative (const Node* var) const { return arg->derivative (var); }
  void optimize (Node** parent)
  { 
    arg -> optimize (&arg); 
    Node* tmp = arg;
    delete (*parent);
    (*parent) = tmp;
    return;
  }
  size_t stackCount(size_t sp, size_t& max_stack) 
  {
    return arg -> stackCount(sp, max_stack);
  }
  // specific
  void addArgument (size_t pos, Node* arg_) { arg = arg_; }
private:
  Node* arg;
};

class NodeMinus : public Node
{
public:
  NodeMinus () : Node(TokenType::UnaryMinus), arg(nullptr) {}
  void deleteTree () { arg->deleteTree (); delete arg; }
  Node* copy () const
  {
    NodeMinus* cp = new NodeMinus();
    cp->arg = arg->copy();
    return cp;
  }
  size_t children () const { return 1; }
  const Node* getChild (size_t n) const { return arg; }
  Node* getChild (size_t n) { return arg; }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const 
  { arg -> evaluate(stack, sp, var, par); stack[sp] *= -1.0; /*std::cout << "evMinus\n";*/ return sp; }
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
  {
    arg -> evaluate(stack, sp, var, par, len);
    for (size_t p = 0; p < len; p++) stack[sp].data[p*stack[sp].skip] *= -1.0;
    return sp;
  }
  void print (std::ostream& out) const 
  {
    P_ERROR_X1( arg != nullptr, "Incerdibly wrong!");
    out << "(";
    out << "-";
    arg->print (out);
    out << ")";
  }
  bool operator != (const Node& node) const 
  { 
    const NodeMinus* nptr = dynamic_cast<const NodeMinus*>(&node);
    if (nptr) return (*arg != *(nptr->arg));
    else return true;
  }
  bool operator < (const Node& node) const 
  {
    if (this->type < node.type) return true; 
    else if (this->type == node.type)
    {
      if (*arg < *(static_cast<const NodeMinus* >(&node)->arg)) return true;
      else return false;
    } 
    return false;
  }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent)
  {
    arg->replaceSymbol (sym, node, &arg);
  }
  void replaceSymbol (const Node& node, Node** parent)
  {
    arg->replaceSymbol (node, &arg);
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
  {
    arg->findFunction (lst, nm);
  }
  Node* derivative (const Node* var) const 
  { 
    NodeMinus* root = new NodeMinus();
    root->arg = arg->derivative (var);
    return root;     
  }
  void optimize (Node** parent)
  { 
    arg -> optimize (&arg); 
    NodeNumber* num = dynamic_cast<NodeNumber*>(arg);
    if (num)
    {
      double val = num->getValue ();
      (*parent)->deleteTree ();
      delete (*parent);
      (*parent) = new NodeNumber(-val);
      return;      
    }
  }
  size_t stackCount(size_t sp, size_t& max_stack) 
  {
    return arg -> stackCount(sp, max_stack);
  }
  // specific
  void addArgument (size_t pos, Node* arg_) { arg = arg_; }
private:
  Node* arg;
};

class NodeEquals : public Node
{
public:
  NodeEquals () : Node(TokenType::Equals), args(2) {}
  void deleteTree () { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; } }
  Node* copy () const
  {
    NodeEquals* cp = new NodeEquals();
    for(size_t k=0; k<args.size(); k++)
    {
      cp->args[k] = args[k]->copy();
    }
    return cp;
  }
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const { return args[n]; }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const 
  {
    P_ERROR_X1( args[0] != nullptr, "Incerdibly wrong!");
    P_ERROR_X1( args[1] != nullptr, "Incerdibly wrong!");
    args[0]->print (out); out << "=";
    args[1]->print (out);
  }
  bool operator != (const Node& node) const 
  { 
    const NodeEquals* nptr = dynamic_cast<const NodeEquals*>(&node);
    if (nptr)
    {
      if (args.size() != nptr->args.size()) return true;
      for(size_t k=0; k<args.size(); k++) if (*args[k] != *(nptr->args[k])) return true;
    }
    else return true;
    return false;
  }
  bool operator < (const Node& node) const 
  { P_MESSAGE1("NOT IMPLEMENTED!\n"); return true; }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent)
  {
    for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
  }
  void replaceSymbol (const Node& node, Node** parent)
  {
    for (size_t k=0; k<args.size(); k++) 
    {
      P_ERROR_X1 (args[k] != nullptr, "terribly wrong");
      args[k]->replaceSymbol (node, &args[k]);
    }
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
  {
    for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
  }
  Node* derivative (const Node* var) const { return new NodeNumber(0.0); }
  void optimize (Node** parent) { for(size_t k=0; k<args.size(); k++) args[k] -> optimize (&args[k]); }  
  size_t stackCount(size_t sp, size_t& max_stack) 
  {
    P_MESSAGE1("NOT IMPLEMENTED!\n"); 
    return sp;
  }
  // specific
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
private:
  std::vector<Node*> args;
};

class NodeSemicolon : public Node
{
public:
  NodeSemicolon () : Node(TokenType::Semicolon) {}
  void deleteTree () { for(auto k=args.begin(); k!=args.end(); k++) { (*k)->deleteTree (); delete *k; } }
  Node* copy () const
  {
    NodeSemicolon* cp = new NodeSemicolon();
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      cp->args.push_back((*k)->copy());
    }
    return cp;
  }
  size_t children () const { return args.size(); }
  const Node* getChild (size_t n) const 
  {
    size_t j=0;
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      if (j == n) return *k;
      j++;
    }
    return nullptr;
  }
  size_t evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  size_t evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
  { P_MESSAGE1("NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const 
  {
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      P_ERROR_X1( (*k) != nullptr, "Incerdibly wrong!");
      (*k)->print (out);
      out << ";\n";
    }
  }
  bool operator != (const Node& node) const { P_MESSAGE1("NOT IMPLEMENTED!\n"); return true; }
  bool operator < (const Node& node) const { P_MESSAGE1("NOT IMPLEMENTED!\n"); return true; }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent)
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->replaceSymbol (sym, node, &(*k));
  }
  void replaceSymbol (const Node& node, Node** parent)
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->replaceSymbol (node, &(*k));
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->findFunction (lst, nm);
  }
  Node* derivative (const Node* var) const { return new NodeNumber(0.0); }
  // specific
  void optimize (Node** parent);
  size_t stackCount(size_t sp, size_t& max_stack) 
  {
    P_MESSAGE1("NOT IMPLEMENTED!\n"); 
    return sp;
  }
  void addArgumentBack (Node* arg) { args.push_back (arg); }
  void addArgumentFront (Node* arg) { args.push_front (arg); }
  void toList(std::list<Node*>& lst)
  {
    lst = std::move(args);
  }
private:
  std::list<Node*> args;
};

void NodeFunction::specialize (NodeFunction** parent)
{
  if (args.size() == 1)
  {
    NodeFunctionA1* fun = new NodeFunctionA1(name);
    fun -> addArgument (0, args[0]);
    delete (*parent);
    (*parent) = fun;
    return;
  }  
  if (args.size() == 2)
  {
    NodeFunctionA2* fun = new NodeFunctionA2(name);
    fun -> addArgument (0, args[0]);
    fun -> addArgument (1, args[1]);
    delete (*parent);
    (*parent) = fun;
    return;
  }
}

void NodeFunction::optimize (Node** parent)
{
  bool isnum = true;
  for (size_t k=0; k<args.size(); k++) 
  {
    args[k]->optimize (&args[k]);
    if (args[k]->type != TokenType::Number) isnum = false;
  }
  // TODO
  // if all are numbers -> evaluate
  if (isnum)
  {
//       int res;
//       static_cast<Node*>(this)->evaluate (&res, std::vector<const int*>(1), std::vector<int>(1));
  }
  // replace power with NodePower
  if (name == "pow")
  {
    NodePower* pw = new NodePower;
    pw -> addArgument (0, args[0]);
    pw -> addArgument (1, args[1]);
    delete (*parent);
    (*parent) = pw;
    return;
  }
  if (*parent == this)
  {
    NodeFunction* tmp = this;
    specialize (&tmp);
    *parent = tmp;
  }
}

void NodeAdd::optimize (Node** parent)
{
  double sum = 0;   // the constant part
  std::vector<Node*> newargs;
  std::vector<double> newmul;
  for (size_t k = 0; k < args.size(); k++)
  {
    args[k]->optimize (&args[k]);
    if (mul[k] == 0.0)
    { 
      args[k]->deleteTree ();
      delete args[k];
      continue;
    }
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[k]);
    if (num)
    {
      sum += mul[k] * (num -> getValue());
      args[k]->deleteTree ();
      delete args[k];
      continue;
    }
    NodeAdd* child = dynamic_cast<NodeAdd*>(args[k]);
    if (child)
    {
      for (size_t p = 0; p < child->children(); p++)
      {
        newargs.push_back (child->args[p]);
        newmul.push_back (mul[k]*child->mul[p]);
      }
      delete child;
      continue;
    }
    NodeMinus* neg = dynamic_cast<NodeMinus*>(args[k]);
    if (neg)
    {
      newargs.push_back (neg->getChild (0));
      newmul.push_back (-1.0*mul[k]);
      delete neg;
      continue;
    }
    NodeTimes* times = dynamic_cast<NodeTimes*>(args[k]);
    if (times)
    {
      newargs.push_back (args[k]);
      newmul.push_back (mul[k]*(times->getMul()));
      times->setMul(1.0);
      continue;
    }
    newargs.push_back (args[k]);
    newmul.push_back (mul[k]);
  }
  args = newargs;
  mul = newmul;
  if (sum != 0.0)
  {
    args.push_back (new NodeNumber(sum));
    mul.push_back (1.0);
  }
  if (args.size() == 0)
  {
    Node* tmp = new NodeNumber (0.0);
    delete *parent;
    *parent = tmp;
    return;
  }
  if (args.size() == 1)
  {
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[0]);
    if (num)
    {
      Node* tmp = new NodeNumber (mul[0] * (num ->getValue()));
      delete num;
      delete *parent;
      *parent = tmp;
      return;
    } else
    {
      if (mul[0] == 1.0)
      {
        Node* tmp = args[0];
        delete *parent;
        *parent = tmp;
        return;
      } else
      {
        // convert into times
        NodeTimes* tmp = new NodeTimes (1);
        tmp -> addArgument (0, args[0], mul[0], false); 
        delete *parent;
        *parent = tmp;
        return;
      }
    }
  }
}

void NodeTimes::optimize (Node** parent)
{
  std::vector<Node*> newargs;
  std::vector<bool> newdivide;
  double newsmul = smul;
  for (size_t k = 0; k < args.size(); k++)
  {
    args[k]->optimize (&args[k]);
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[k]);
    if (num)
    {
      if (divide[k]) newsmul /= num -> getValue();
      else newsmul *= num -> getValue();
      args[k]->deleteTree ();
      delete args[k];
      continue;
    }
    NodeTimes* child = dynamic_cast<NodeTimes*>(args[k]);
    if (child)
    {
      args[k] = child->args[0];
      if (divide[k]) newsmul /= child->smul; 
      else newsmul *= child->smul;
      for (size_t p = 0; p < child->children(); p++)
      {
        newargs.push_back (child->args[p]);
        if (divide[k]) newdivide.push_back (!child->divide[p]);
        else newdivide.push_back (child->divide[p]);
      }
      delete child;
      continue;
    }
    NodeMinus* neg = dynamic_cast<NodeMinus*>(args[k]);
    if (neg)
    {
      newsmul *= -1.0;
      newargs.push_back (neg->getChild (0));
      newdivide.push_back (child->divide[k]);
      delete neg;
      continue;
    }
    newargs.push_back (args[k]);
    newdivide.push_back (divide[k]);
  }
  args = newargs;
  divide = newdivide;
  smul = newsmul;
  if (smul == 0.0)
  {
    (*parent)->deleteTree ();
    delete *parent;
    *parent = new NodeNumber(0.0);
    return;
  }
  // single argument
  if ((smul == 1.0)&&(args.size() == 1)&&(!divide[0]))
  {
    // self destruct
    Node* tmp = args[0];
    delete *parent;
    *parent = tmp;
    return;
  }
  if (args.size() == 0)
  {
    const double smul_copy = smul;
    delete *parent;
    *parent = new NodeNumber(smul_copy);
    return;    
  }
}

void NodeSemicolon::optimize (Node** parent)
{
  auto it=args.begin();
  while (it != args.end())
  {
    (*it)->optimize (&(*it));
    NodeSemicolon* child = dynamic_cast<NodeSemicolon*>(*it);
    if (child)
    {
      for (auto p=child->args.begin(); p != child->args.end(); p++)
      {
        (*p)->optimize (&(*p));
        args.insert (it, *p);
      }
      it = args.erase (it);
      delete child;
    } else
    it++;
  }
}

// DERIVATIVE IMPLEMENTATIONS

Node* NodeNumber::derivative (const Node* var) const 
{
  return new NodeNumber(0.0);
}

Node* NodeSymbol::derivative (const Node* var) const
{
  if (var->type == TokenType::Symbol)
  {
    const NodeSymbol* node = static_cast<const NodeSymbol* >(var);
    {
      if (node->name == name)
      {
        return new NodeNumber(1.0);
      }
    }
  }
  return new NodeNumber(0.0);
}

Node* NodeVar::derivative (const Node* var) const
{
  if (var->type == TokenType::Var)
  {
    const NodeVar* node = static_cast<const NodeVar* >(var);
    {
      if (node->idx == idx)
      {
        return new NodeNumber(1.0);
      }
    }
  }
  return new NodeNumber(0.0);
}

Node* NodePar::derivative (const Node* var) const
{
  if (var->type == TokenType::Par)
  {
    const NodePar* node = static_cast<const NodePar* >(var);
    {
      if (node->idx == idx)
      {
        return new NodeNumber(1.0);
      }
    }
  }
  return new NodeNumber(0.0);
}

Node* NodeAdd::derivative (const Node* var) const
{
  NodeAdd* cp = new NodeAdd(args.size());
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->derivative (var);
    cp->mul[k] = mul[k];
  }
  return cp;
}

Node* NodeTimes::derivative (const Node* var) const  
{
  // f = (f1 f2 f2 ..) / (g1 g2 g3 ..)
  NodeAdd* root = new NodeAdd(args.size());
  for(size_t k=0; k<args.size(); k++)
  {
    NodeTimes* cp = new NodeTimes(0);
    cp->smul = smul;
    if (divide[k])
    {
      // divisions
      cp->smul *= -1.0;
      cp->args.push_back (args[k]->derivative (var));
      cp->divide.push_back (false);
      cp->args.push_back (args[k]->copy ());
      cp->divide.push_back (true);
      cp->args.push_back (args[k]->copy ());
      cp->divide.push_back (true);
    }
    else
    {
      cp->args.push_back (args[k]->derivative (var));
      cp->divide.push_back (false);
    }
    for (size_t l=0; l<args.size(); l++)
    { 
      // not the ones being differentiated
      if (k != l) 
      { 
        cp->args.push_back (args[l]->copy());
        cp->divide.push_back (divide[l]);
      }
    }
    root -> addArgument (k, cp, 1.0);
  }
  return root;
}

Node* NodeFunction::derivative (const Node* var) const
{
  NodeTimes* root = new NodeTimes(2);
  if (name == "abs")
  {
    // -> sign
    NodeFunctionA1* fun = new NodeFunctionA1("sign");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
  }
  if (name == "exp")
  {
    // -> exp
    NodeFunctionA1* fun = new NodeFunctionA1("exp");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);    
  }
  if (name == "log") 
  {
    // -> 1/x
    root -> addArgument (0, args[0] -> copy (), 1.0, true);
  }
  if (name == "sin")
  {
    // -> cos
    NodeFunctionA1* fun = new NodeFunctionA1("cos");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);    
  }
  if (name == "cos")
  {
    // -> -sin
    NodeFunctionA1* fun = new NodeFunctionA1("sin");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, -1.0, false);    
  }
  if (name == "tan")
  {
    // -> cos^(-2)
    NodePower* pow = new NodePower;
    NodeFunctionA1* fun = new NodeFunctionA1("cos");
    fun -> addArgument (0, args[0] -> copy ());
    pow -> addArgument (0, fun);
    pow -> addArgument (1, new NodeNumber(2.0));
    root -> addArgument (0, pow, 1.0, true);    
  }
  if (name == "asin")
  {
    // -> 1/{(1-x^2)^(1/2)}
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (0, new NodeNumber(1.0), 1.0);
    add -> addArgument (1, pow_square, -1.0);   
    NodePower* pow_sqrt = new NodePower;
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber(0.5));
    root -> addArgument (0, pow_sqrt, 1.0, true);
  }
  if (name == "acos")
  {
    // -> -1/{(1-x^2)^(1/2)}
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (0, new NodeNumber(1.0), 1.0);
    add -> addArgument (1, pow_square, -1.0);   
    NodePower* pow_sqrt = new NodePower;
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber(0.5));
    root -> addArgument (0, pow_sqrt, -1.0, true);
  }
  if (name == "atan")
  {
    // -> 1/(1-x^2)
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (0, new NodeNumber(1.0), 1.0);
    add -> addArgument (1, pow_square, -1.0);   
    root -> addArgument (0, add, 1.0, true);
  }
  if (name == "sinh")
  {
    // -> cosh
    NodeFunctionA1* fun = new NodeFunctionA1("cosh");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);    
  }
  if (name == "cosh")
  {
    // -> sinh
    NodeFunctionA1* fun = new NodeFunctionA1("sinh");
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);    
  }
  if (name == "tanh")
  {
    // -> cosh^(-2)
    NodePower* pow = new NodePower;
    NodeFunctionA1* fun = new NodeFunctionA1("cosh");
    fun -> addArgument (0, args[0] -> copy ());
    pow -> addArgument (0, fun);
    pow -> addArgument (1, new NodeNumber(2.0));
    root -> addArgument (0, pow, 1.0, true);    
  }
  if (name == "asinh")
  {
    // -> 1/(1+x^2)^(1/2)
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (0, new NodeNumber(1.0), 1.0);
    add -> addArgument (1, pow_square, 1.0);   
    NodePower* pow_sqrt = new NodePower;
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber(0.5));
    root -> addArgument (0, pow_sqrt, 1.0, true);
  }
  if (name == "acosh")
  {
    // -> 1/(x^2-1)^(1/2)
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (1, new NodeNumber(1.0), -1.0);
    add -> addArgument (0, pow_square, 1.0);   
    NodePower* pow_sqrt = new NodePower;
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber(0.5));
    root -> addArgument (0, pow_sqrt, 1.0, true);
  }
  if (name == "atanh")
  {
    // -> 1/(1-x^2)
    NodePower* pow_square = new NodePower;
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber(2.0));
    NodeAdd* add = new NodeAdd(2);
    add -> addArgument (0, new NodeNumber(1.0), 1.0);
    add -> addArgument (1, pow_square, -1.0);   
    root -> addArgument (0, add, 1.0, true);
  }
  if (name == "pow")
  {
    // promote to NodePower
    delete root;
    NodePower* pw = new NodePower;
    pw -> addArgument (0, args[0]);
    pw -> addArgument (1, args[1]);
    Node* retval = pw -> derivative (var);
    delete pw;
    return retval;
  }
  if (name == "sign")
  {
    delete root;
    return new NodeNumber(0.0);
  }
  if (name == "heaviside")
  {
    delete root;
    return new NodeNumber(0.0);
  }
  root -> addArgument (1, args[0] -> derivative (var), 1.0, false);
  return root;
}

// Evaluating the expression

size_t NodeNumber::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  stack[sp] = value;
  return sp;
}

size_t NodeSymbol::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  P_MESSAGE2("Undefined symbol : ", name);
  stack[sp] = 0;
  return sp;
}

size_t NodeVar::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  stack[sp] = var[idx];
  return sp;
}

size_t NodePar::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  stack[sp] = par[idx];
  return sp;
}

size_t NodeFunction::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  P_MESSAGE3("NodeFunction::evaluate : it is a base class, cannot be called ", name, "\n");
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k, var, par); }
  if (spi != sp + args.size() - 1) 
    P_MESSAGE7("NodeFunction::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size()-1, "\n");
  stack[sp] = 1.0;
  return sp;
}

size_t NodeFunctionA1::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  size_t sp1 = args[0] -> evaluate(stack, sp, var, par);
  stack[sp1] = function(stack[sp1]);
  return sp1;
}

size_t NodeFunctionA2::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  size_t sp1 = args[0] -> evaluate(stack, sp, var, par);
  size_t sp2 = args[1] -> evaluate(stack, sp+1, var, par);
  stack[sp1] = function(stack[sp1], stack[sp2]);
  return sp1;
}

size_t NodeAdd::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, var, par); }
  if (spi != sp + args.size()) P_MESSAGE7("NodeAdd::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
  stack[sp] = 0;
  for (size_t k=0; k<args.size(); k++)
  {
    stack[sp] += mul[k] * stack[sp+k+1];
  }
  return sp;
}

size_t NodeTimes::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, var, par); }
  if (spi != sp + args.size()) P_MESSAGE7("NodeTimes::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
  stack[sp] = smul;
  for (size_t k=0; k<args.size(); k++)
  {
    if (divide[k]) stack[sp] /= stack[sp+k+1];
    else stack[sp] *= stack[sp+k+1];
  }
  return sp;
}

size_t NodePower::evaluate (std::vector<double>& stack, size_t sp, const std::vector<double>& var, const std::vector<double>& par) const
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k, var, par); }
  if (spi != sp + args.size() - 1) P_MESSAGE7("NodePower::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size()-1, "\n");
  stack[sp] = pow (stack[sp], stack[sp+1]);
  return sp;
}

// ----------------------------------------------------------------------------------------------------------------------------
// Vector arguments

size_t NodeNumber::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = value;
  }
  return sp;
}

size_t NodeSymbol::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  P_MESSAGE2("Undefined symbol : ", name);
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = 0;
  }
  return sp;
}

size_t NodeVar::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = var[idx].data[k*var[idx].skip];
  }
  return sp;
}

size_t NodePar::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = par[idx];
  }
  return sp;
}

size_t NodeFunction::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  P_MESSAGE3("NodeFunction::evaluate : it is a base class, cannot be called ", name, "\n");
  return sp;
}

size_t NodeFunctionA1::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  args[0] -> evaluate(stack, sp, var, par, len);
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = this->function(data[skip*k]);
  }
  return sp;
}

size_t NodeFunctionA2::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  args[0] -> evaluate(stack, sp, var, par, len);
  args[1] -> evaluate(stack, sp+1, var, par, len);
  double* data1 = stack[sp].data;
  size_t skip1 = stack[sp].skip;
  double* data2 = stack[sp+1].data;
  size_t skip2 = stack[sp+1].skip;
  for (size_t k = 0; k < len; k++)
  {
    data1[skip1*k] = this->function(data1[skip1*k], data2[skip2*k]);
  }
  return sp;
}

size_t NodeAdd::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, var, par, len); }
  if (spi != sp + args.size()) P_MESSAGE7("NodeAdd::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
  for (size_t p = 0; p < len; p++)
  {
    stack[sp].data[p*stack[sp].skip] = 0;
  }
  for (size_t k = 0; k < args.size(); k++)
  {
    for (size_t p = 0; p < len; p++)
    {
      stack[sp].data[p*stack[sp].skip] += mul[k] * stack[sp+k+1].data[p*stack[sp+k+1].skip];
    }
  }
  return sp;
}

size_t NodeTimes::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, var, par, len); }
  if (spi != sp + args.size()) P_MESSAGE7("NodeTimes::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
  for (size_t p = 0; p < len; p++)
  {
    stack[sp].data[p*stack[sp].skip] = smul;
  }
  for (size_t k = 0; k < args.size(); k++)
  {
    for (size_t p = 0; p < len; p++)
    {
      if (divide[k]) stack[sp].data[p*stack[sp].skip] /= stack[sp+k+1].data[p*stack[sp+k+1].skip];
      else stack[sp].data[p*stack[sp].skip] *= stack[sp+k+1].data[p*stack[sp+k+1].skip];
    }
  }
  return sp;
}

size_t NodePower::evaluate (std::vector<Value>& stack, size_t sp, const std::vector<Value>& var, const std::vector<double>& par, size_t len) const
{
  size_t spi = sp;
  spi = args[0] -> evaluate(stack, sp, var, par, len);
  spi = args[1] -> evaluate(stack, sp+1, var, par, len);
  if (spi != sp + args.size() - 1) P_MESSAGE7("NodePower::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size()-1, "\n");
  for (size_t p = 0; p < len; p++)
  {
    stack[sp].data[p*stack[sp].skip] = pow (stack[sp].data[p*stack[sp].skip], stack[sp+1].data[p*stack[sp+1].skip]);
  }
  return sp;
}

// using namespace ExpTree;

// There are a number of limits to tokens
// alphanumeric is an identifier
// operators inbetween are the delimiters and they are tokens separately
// 

enum class CharType : int { Op, Identifier, Space, Invalid };

static inline CharType findCharType (char c)
{
  if ( isalnum ( c ) || (c == '_') || (c == '.') ) return CharType::Identifier;
  else if ( (c == '(') || (c == ')') || (c == '*') || (c == '/') || (c == '+') || (c == '-') || (c == '^') || (c == ',') || (c == '=') || (c == ';')) return CharType::Op;
  else if ( isspace ( c ) ) return CharType::Space;
  else return CharType::Invalid;
}

void tokenize (std::vector<Token>& stream, const std::map<TokenType, TokenPrecedence>& token_table, const std::string& expression)
{
  const char *str = expression.c_str ();
  const size_t len = expression.size ();
  
  size_t k = 0;
  while (k < len)
  {
    if (findCharType (str[k]) == CharType::Space)
    {
      k++;
      continue;
    } else if (findCharType (str[k]) == CharType::Identifier)
    {
      if (std::isdigit(str[k]) || str[k] == '.')
      {
        char* endptr;
        strtod (&str[k], &endptr);
        size_t p = (endptr - &str[k])/sizeof(char);
        std::string tn (&str[k], p);
        k = k + p;
        stream.push_back (Token (TokenType::Symbol, tn) );
        continue;
      } else
      {
        size_t p = 0;
        while (findCharType (str[k+p]) == CharType::Identifier) p++;
        std::string tn (&str[k], p);
        k = k + p;
        stream.push_back (Token (TokenType::Symbol, tn) );
        continue;
      }
    } else if (findCharType (str[k]) == CharType::Op)
    {
      std::string tn (&str[k],1);
      // a single character
      switch (tn[0])
      {
        case '+':
          stream.push_back (Token (TokenType::Plus, tn) );
          break;
        case '-':
          stream.push_back (Token (TokenType::Minus, tn) );
          break;
        case '*':
          stream.push_back (Token (TokenType::Times, tn) );
          break;
        case '/':
          stream.push_back (Token (TokenType::Divide, tn) );
          break;
        case '^':
          stream.push_back (Token (TokenType::Power, tn) );
          break;
        case '(':
          stream.push_back (Token (TokenType::LeftParen, tn) );
          break;
        case ')':
          stream.push_back (Token (TokenType::RightParen, tn) );
          break;
        case ',':
          stream.push_back (Token (TokenType::Comma, tn) );
          break;
        case '=':
          stream.push_back (Token (TokenType::Equals, tn) );
          break;
        case ';':
          stream.push_back (Token (TokenType::Semicolon, tn) );
          break;
        default:
          P_MESSAGE3("unknown operator ", tn, "\n");
      }
      k++;
    } else
    {
      P_MESSAGE3("Invalid character ", str[k], "\n");
      k++;
    }
  }
  if (stream[0].type == TokenType::Minus) stream[0].type = TokenType::UnaryMinus;
  if (stream[0].type == TokenType::Plus) stream[0].type = TokenType::UnaryPlus;
  for (size_t k=1; k < stream.size(); ++k)
  {
    // identify functions
    if ((stream[k].type == TokenType::LeftParen)&&(stream[k-1].type == TokenType::Symbol))
    {
      stream[k-1].type = TokenType::Function;
    }
    // identify unary minus
    if ((stream[k].type == TokenType::Minus)&&(stream[k-1].type != TokenType::Symbol)&&(stream[k-1].type != TokenType::RightParen))
    {
      stream[k].type = TokenType::UnaryMinus;
    }
    // identify unary plus
    if ((stream[k].type == TokenType::Plus)&&(stream[k-1].type != TokenType::Symbol)&&(stream[k-1].type != TokenType::RightParen))
    {
      stream[k].type = TokenType::UnaryPlus;
    }
  }
  // check if the last token is a semicolon
  if (stream.back().type == TokenType::Semicolon) stream.pop_back();
  else P_MESSAGE1("Missing semicolon\n");
}

void opToTree (std::list<Node*>& ts, const Token& tk, size_t nargs)
{
  if (ts.size() >= nargs)
  {
    switch (tk.type)
    {
      case TokenType::Symbol:
      {
        char* endptr;
        double res = strtod (tk.name.c_str(), &endptr);
        if (endptr != tk.name.c_str())
        {
          NodeNumber* node = new NodeNumber (res);
          ts.push_back (node);          
        } else
        {
          NodeSymbol* node = new NodeSymbol (tk.name);
          ts.push_back (node);
        }
      }
      break;
      case TokenType::Function:
      {
        if (nargs == 1)
        {
          NodeFunctionA1* node = new NodeFunctionA1 (tk.name);
          // put the argument from the stack
          node->addArgument (0, ts.back());
          ts.pop_back ();
          ts.push_back (node);
        } else if (nargs == 2)
        {
          NodeFunctionA2* node = new NodeFunctionA2 (tk.name);
          // put the argument from the stack
          node->addArgument (1, ts.back());
          ts.pop_back ();
          node->addArgument (0, ts.back());
          ts.pop_back ();
          ts.push_back (node);
        } else
        {
          NodeFunction* node = new NodeFunction (tk.name, nargs);
          // put the argument from the stack
          for (size_t k=0; k<nargs; k++)
          {
             node->addArgument (nargs - k - 1, ts.back());
             ts.pop_back ();
          }
          ts.push_back (node);
        }
      }
      break;
      case TokenType::UnaryPlus:
      {
        // just ignore it
      }
      break;
      case TokenType::UnaryMinus:
      {
        NodeMinus* node = new NodeMinus ();
        node->addArgument (0, ts.back());
        ts.pop_back ();
        ts.push_back (node);
      }
      break;
      case TokenType::Plus:
      {
        NodeAdd* node = new NodeAdd (nargs);
        // put the argument from the stack
        for (size_t k=0; k<nargs; k++)
        {
          node->addArgument (nargs - k - 1, ts.back(), 1.0);
          ts.pop_back ();
        }
        node->sort ();
        ts.push_back (node);
      }
      break;
      case TokenType::Minus:
      {
        NodeAdd* node = new NodeAdd (nargs);
        // put the argument from the stack
        if (nargs == 2)
        {
          node->addArgument (1, ts.back(), -1.0);
          ts.pop_back ();
          node->addArgument (0, ts.back(), 1.0);
          ts.pop_back ();
        } else P_MESSAGE1("Wrong no. args\n");
        node->sort ();
        ts.push_back (node);
      }
      break;
      case TokenType::Times:
      {
        NodeTimes* node = new NodeTimes (nargs);
        // put the argument from the stack
        for (size_t k=0; k<nargs; k++)
        {
          node->addArgument (nargs - k - 1, ts.back(), 1.0, false);
          ts.pop_back ();
        }
        node->sort ();
        ts.push_back (node);
      }
      break;
      case TokenType::Divide:
      {
        NodeTimes* node = new NodeTimes (nargs);
        // put the argument from the stack
        for (size_t k=0; k<nargs; k++)
        {
          if (k == 0) node->addArgument (nargs - k - 1, ts.back(), 1.0, true);
          else node->addArgument (nargs - k - 1, ts.back(), 1.0, false);
          ts.pop_back ();
        }
        node->sort ();
        ts.push_back (node);
      }
      break;
      case TokenType::Power:
      {
        // nargs == 2 ??
        NodePower* node = new NodePower ();
        // put the argument from the stack
        node->addArgument (1, ts.back());
        ts.pop_back ();
        node->addArgument (0, ts.back());
        ts.pop_back ();
        ts.push_back (node);
      }
      break;
      case TokenType::Equals:
      {
        NodeEquals* node = new NodeEquals ();
        // put the argument from the stack
        node->addArgument (1, ts.back());
        ts.pop_back ();
        node->addArgument (0, ts.back());
        ts.pop_back ();
        ts.push_back (node);
      }
      break;
      case TokenType::Semicolon:
      {
        NodeSemicolon* node = new NodeSemicolon ();
        // put the argument from the stack
        Node* nd2 = ts.back();
        ts.pop_back ();
        Node* nd1 = ts.back();
        ts.pop_back ();
        NodeSemicolon* semi = dynamic_cast<NodeSemicolon*>(nd1);
        if (semi != nullptr)
        {
          for (size_t k = 0; k < semi -> children(); k++)
          {
            node->addArgumentBack (const_cast<Node*>(semi->getChild(k)));
          }
          delete semi;
        } else node->addArgumentBack (nd1);
        // second argument
        semi = dynamic_cast<NodeSemicolon*>(nd2);
        if (semi != nullptr)
        {
          for (size_t k = 0; k < semi -> children(); k++)
          {
            node->addArgumentBack (const_cast<Node*>(semi->getChild(k)));
          }
          delete semi;
        } else node->addArgumentBack (nd2);
        ts.push_back (node);
      }
      break;
      default:
        P_MESSAGE1("Unexpected operator\n");
        break;
    }
  } else P_MESSAGE1("Too short stack\n");
}

Node* toPostfix (const std::vector<Token>& token_stream, const std::map<TokenType, TokenPrecedence>& token_table)
{
  std::list<Token> ops;
  std::list<size_t> funargs;
  std::list<Node*> treestack;
  size_t k = 0;
  while (k < token_stream.size())
  {
    if (token_stream[k].type == TokenType::Symbol)
    {
      if (!funargs.empty()) if (funargs.back() == 0) funargs.back() = 1;
      // bind items on the stack ...
      opToTree (treestack, token_stream[k], 0);
    }
    else if (token_stream[k].type == TokenType::Function) 
    {
      ops.push_back (token_stream[k]);
    }
    else if (token_stream[k].type == TokenType::Comma)
    {
      if (!funargs.empty()) funargs.back() += 1;
      else P_MESSAGE1("not a function\n");
      while ( (!ops.empty())&&(ops.back().type != TokenType::LeftParen) )
      {
        // bind items on the stack ...
        opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
        ops.pop_back();
      }
      if (ops.empty()) P_MESSAGE1("1Misplaced comma or missing \')\'\n"); 
    }
    else if (token_stream[k].isOp())
    {
      while (!ops.empty())
      {
        if ( (ops.back().isOp())&&
             (((token_table.at(token_stream[k].type).leftAssoc)&&
              (token_table.at(token_stream[k].type).precedence >= token_table.at(ops.back().type).precedence))||
              (token_table.at(token_stream[k].type).precedence > token_table.at(ops.back().type).precedence)) )
        {
          // bind items on the stack ...
          opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
          ops.pop_back();
        } else break;
      }
      ops.push_back(token_stream[k]);
    }
    else if (token_stream[k].type == TokenType::LeftParen)
    {
      funargs.push_back(0);
      ops.push_back(token_stream[k]);
    }
    else if (token_stream[k].type == TokenType::RightParen)
    {
      while (!ops.empty())
      {
        if (ops.back().type != TokenType::LeftParen)
        {
          // bind items on the stack ...
          opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
          ops.pop_back();
        } else 
        { 
          break;
        }
      } 
      // the right paren is still on the stack
      if (!ops.empty()) ops.pop_back();
      else
      {
        // there was noting on the stack
	std::ostringstream msg;
	msg << "2Mismatched prenthesis\n BEGIN@ ";
        for (size_t p=0; p<k+1; p++) msg << token_stream[p].name;
        msg << " @ ";
        for (size_t p=k+1; p<token_stream.size(); p++) msg << token_stream[p].name;
        msg << " @END\n ";
        P_MESSAGE1(msg.str().c_str());
      }
      if (!ops.empty()) if (ops.back().type == TokenType::Function)
      {
        // bind items on the stack ...
        opToTree (treestack, ops.back(), funargs.back());
        ops.pop_back();
      }
      funargs.pop_back();
    }
    k++;
  }
  // no more tokens
  while (!ops.empty())
  {
    if (ops.back().type == TokenType::LeftParen) P_MESSAGE1("Mismatched prenthesis\n");
    else
    {
      // bind items on the stack ...
      opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
      ops.pop_back();
    }
  }
  // finished
  if (treestack.size() != 1)
  {
    std::ostringstream msg;
    msg << "Operations remaining " << treestack.size() << "\n";
    for (auto k=treestack.begin(); k != treestack.end(); ++k)
    {
      msg << "@@OP =\n\t";
      if (*k) (*k)->print (msg); msg << "\n";
    }
    P_MESSAGE1(msg.str().c_str());
  }
  return treestack.front ();
}

// needs to handle period(), time(), PAR_PERIOD
void splitExpression (std::string& sysName, 
  std::vector<NodeSymbol*>& var_name,
  std::vector<Node*>& var_dot,
  std::vector<Node*>& var_init, 
  std::vector<NodeSymbol*>& par_name,
  std::vector<Node*>& par_value,
  std::vector<NodeSymbol*>& time,
  std::vector<Node*>& macro,
  Node* root)
{
  // set this up for the period
  par_name.push_back(nullptr);
  par_value.push_back(nullptr);
  if (root->type == TokenType::Semicolon)
  {
    Node* cp = root -> copy ();
    std::list<Node*> lst;
    static_cast<NodeSemicolon*>(cp)->toList(lst);
    for (auto it : lst)
    {
      if (it->type == TokenType::Equals)
      {
        const Node* left = it->getChild (0);
        const Node* right = it->getChild (1);
        
        if (left->type == TokenType::Function)
        {
          const NodeFunction* fun = static_cast<const NodeFunction*>(left);
          if (fun->getName() == "dot")
          {
            // state variable
            if (fun->children() == 1)
            {
              const Node* statevariable = fun->getChild(0);
              if (statevariable->type == TokenType::Symbol)
              {
                var_name.push_back (static_cast<NodeSymbol*>(statevariable->copy()));
                var_dot.push_back (right->copy());
              }
            } else P_MESSAGE1("Dot: wrong number of arguments\n");
          } else if (fun->getName() == "par")
          {
            // parameter
            if (fun->children() == 1)
            {
              const Node* par = fun->getChild(0);
              if (par->type == TokenType::Symbol)
              {
                par_name.push_back (static_cast<NodeSymbol*>(par->copy()));
                par_value.push_back (right->copy());
              }
            } else P_MESSAGE1("Dot: wrong number of arguments\n");
          } else if (fun->getName() == "init")
          {
            // starting solution
            if (fun->children() == 1)
            {
              const Node* init = fun->getChild(0);
              if (init->type == TokenType::Symbol)
              {
                const NodeSymbol* sym = static_cast<const NodeSymbol*>(init);
                // find the corresponding state variable
                var_init.resize(var_name.size());
                for (size_t k=0; k<var_name.size(); k++)
                {
                  if (var_name[k]->getName() == sym->getName())
                  {
                    var_init[k] = right->copy();
                    break;
                  }
                }
              }
            } else P_MESSAGE1("Dot: wrong number of arguments\n");
          } else if (fun->getName() == "period")
          {
            if (par_name.size() > 0) P_ERROR_X1(par_name[0] == nullptr, "Period is already defined.");
            P_ERROR_X1(fun->children() == 1, "Dot: wrong number of arguments");
            P_ERROR_X1(fun->getChild(0)->type == TokenType::Symbol, "Not a symbol");
            par_name[0] = static_cast<NodeSymbol*>(fun->getChild(0)->copy());
            par_value[0] = right->copy();
          } else if (fun->getName() == "time")
          {
            P_ERROR_X1(fun->children() == 0, "Too many arguments");
            P_ERROR_X1(right->type == TokenType::Symbol, "Not a symbol");
            time.push_back (static_cast<NodeSymbol*>(right->copy()));
          } else if (fun->getName() == "vfname")
          {
            P_ERROR_X1(fun->children() == 0, "Too many arguments");
            P_ERROR_X1(right->type == TokenType::Symbol, "Not a symbol");
            sysName = static_cast<const NodeSymbol*>(right)->getName ();
          } else
          {
            // a simple macro
            macro.push_back(it->copy());
          }
        } else if (left->type == TokenType::Symbol)
        {
          macro.push_back(it->copy());
        } else 
        {
          left->print (std::cout); std::cout << " == "; right->print (std::cout); std::cout << "\n";
          std::cout.flush();
          P_MESSAGE1("Not a definition\n");
        }
      } else
      {
        it->print (std::cout); std::cout << "\n";
        P_MESSAGE1("Not a definition\n");
      }
    }
    for (auto it : lst) { it->deleteTree (); delete it; }
    cp -> deleteTree ();
    delete cp;
  } else
  {
    P_MESSAGE1("Cannot split expression\n");
  }
  if (par_name[0] == nullptr)
  {
    par_name[0] = new NodeSymbol("T");
    par_value[0] = new NodeNumber(1.0);
  }
}

class ExpressionParser
{
public:
  ExpressionParser ()
  {
    token_table[TokenType::LeftParen] = { 1, true, 0 };
    token_table[TokenType::RightParen] = { 1, true, 0 };
    token_table[TokenType::Function] = { 2, true, 0 };
    token_table[TokenType::UnaryPlus] = { 3, true, 1 };
    token_table[TokenType::UnaryMinus] = { 3, true, 1 };
    token_table[TokenType::Power] = { 4, false, 2 };
    token_table[TokenType::Times] = { 5, true, 2 };
    token_table[TokenType::Divide] = { 5, true, 2 };
    token_table[TokenType::Plus] = { 7, true, 2 };
    token_table[TokenType::Minus] = { 7, true, 2 };
    token_table[TokenType::Comma] = { 9, true, 0 };
    token_table[TokenType::Equals] = { 10, false, 2 };
    token_table[TokenType::Semicolon] = { 11, true, 2 };
  }
  ~ExpressionParser () {}
  Node* parse (const std::string& str)
  {
    std::vector<Token> token_stream;
    tokenize (token_stream, token_table, str);
    return toPostfix (token_stream, token_table);
  }
private:
  std::map<TokenType, TokenPrecedence> token_table;
};

Expression::Expression (const std::string& str)
{
  ExpressionParser parser;
  root = parser.parse (str);
}

void Expression::fromString (const std::string& str)
{
  if (root != nullptr)
  {
    root -> deleteTree ();
    delete root;
    root = nullptr;
  }
  ExpressionParser parser;
  root = parser.parse (str);
}

void Expression::knutSplit (
  std::string& sysName,
  std::vector<std::string>& varName,
  std::vector<Expression>& varDotExpr,
  std::vector<Expression>& varInit,
  std::vector<Expression>& delayExpr,
  std::vector<std::string>& parName,
  std::vector<double>& parInit )
{
  std::vector<NodeSymbol*> var_name;
  std::vector<Node*> var_dot;
  std::vector<Node*> var_init; 
  std::vector<NodeSymbol*> par_name;
  std::vector<Node*> par_value;
  std::vector<NodeSymbol*> time;
  std::vector<Node*> macro;

  // Adding Pi
  NodeEquals* pi_eq = new NodeEquals;
  pi_eq -> addArgument (0, new NodeSymbol("Pi"));
  pi_eq -> addArgument (1, new NodeNumber(M_PI));
  macro.push_back (pi_eq);

  splitExpression (sysName, var_name, var_dot, var_init, par_name, par_value, time, macro, root);

  // if no independent variable or period is set we set the default
  if (time.empty()) time.push_back (new NodeSymbol("t"));
  NodeVar* time_var = new NodeVar(0);

  // adding internal parameters
  par_name.push_back (new NodeSymbol ("PAR_ANGLE"));
  par_value.push_back (new NodeNumber (0.0));
  par_name.push_back (new NodeSymbol ("PAR_PERIOD"));
  par_value.push_back (new NodeNumber (1.0));
  par_name.push_back (new NodeSymbol ("PAR_ROT_NUM"));
  par_value.push_back (new NodeNumber (0.0));
    
  std::vector<NodeVar*> var_idx;
  std::vector<NodePar*> par_idx;
  
  // dummy variables and parameters
  std::vector<double> var_numval;

  varName.resize (var_name.size());
  for (size_t k = 0; k < var_name.size(); k++)
  {
    var_idx.push_back (new NodeVar(1 + k)); // idx == 0 is the time
    varName[k] = var_name[k]->getName ();
  }

  parName.resize (par_name.size());
  for (size_t k = 0; k < par_name.size(); k++)
  {
    par_idx.push_back (new NodePar(k)); // idx == 0 is the period
    parName[k] = par_name[k]->getName ();
  }

  // Applying the macros to both parameters and variables
  for (size_t p = 0; p < macro.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      if (p != q) macro[p] -> replaceSymbol (*macro[q], &macro[p]);
    }
  }

  for (size_t p = 0; p < var_name.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      var_dot[p] -> replaceSymbol (*macro[q], &var_dot[p]);
    }
  }
  
  for (size_t p = 0; p < var_init.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      var_init[p] -> replaceSymbol (*macro[q], &var_init[p]);
    }
    for (size_t q = 0; q < par_name.size(); q++)
    {
      var_init[p] -> replaceSymbol (*par_name[q], *par_value[q], &var_init[p]);
    }
    var_init[p] -> replaceSymbol (*time[0], *time_var, &var_init[p]);
  }

  for (size_t p = 0; p < par_name.size(); p++)
  {
    // other params
    for (size_t q = 0; q < par_name.size(); q++)
    {
      if (p != q) par_value[p] -> replaceSymbol (*par_name[q], *par_value[q], &par_value[p]);
    }
    // twice
    for (size_t q = 0; q < par_name.size(); q++)
    {
      if (p != q) par_value[p] -> replaceSymbol (*par_name[q], *par_value[q], &par_value[p]);
    }
    for (size_t q = 0; q < macro.size(); q++)
    {
      par_value[p] -> replaceSymbol (*macro[q], &par_value[p]);
    }
  }

  // fill up the parameter values
  parInit.resize (par_name.size());
  for (size_t p = 0; p < par_name.size(); p++)
  {
    size_t stmax = 0;
    par_value[p] -> stackCount(0,stmax);
    std::vector<double> stack(stmax);
    par_value[p] -> evaluate(stack, 0, var_numval, parInit);
    parInit[p] = stack[0];
  }

  // fill up the starting point
  varInit.resize (var_init.size());
  for (size_t p = 0; p < var_init.size(); p++)
  {
    varInit[p].fromNode (var_init[p]);
  }
  
  // START :: HANDLING THE DELAYS
  // do not delete any of these since they are owned by the tree!
  std::vector<const NodeFunction*> delays;
  for (size_t p = 0; p < var_name.size(); p++)
  {
    var_dot[p] -> findFunction (delays, "delay");
  }
  // sort the delays
  std::sort (delays.begin(), delays.end(), 
    [&](const NodeFunction* a, const NodeFunction* b) { return (*a->getChild (1)) < (*b->getChild (1)); } );

  // keep the delays
  std::vector<const NodeFunction*> unique_delays;
  std::vector<size_t> delay_id(delays.size());
  
  // remove the extra delays
  if (!delays.empty()) 
  {
    unique_delays.push_back (delays[0]);
    delay_id[0] = 0;
  }
  size_t idx = 0;
  for (size_t k = 1; k < delays.size(); k++)
  {
    if (*(delays[k]->getChild (1)) != *(unique_delays[idx]->getChild (1)))
    {
      idx++;
      unique_delays.push_back (delays[k]);
    }
    delay_id[k] = idx;
  }

  // transforming vars into delayed values
  for (size_t p = 0; p < delays.size(); p++)
  {
    std::vector<NodeVar*> var_del(var_name.size());
    NodeFunction* fun = const_cast<NodeFunction*>(delays[p]);
    for (size_t q = 0; q < var_name.size(); q++)
    {
      NodeVar var_del(1 + q + var_name.size()*(delay_id[p]+1));
      fun->replaceSymbol (*var_name[q], var_del, nullptr);
    }
  }
  
  // fill up the delays
  delayExpr.resize (unique_delays.size());
  for (size_t r = 0; r < unique_delays.size(); r++)
  {
    Node* delay = unique_delays[r]->getChild (1)->copy();
    // replace parameters
    for (size_t q = 0; q < par_name.size(); q++)
    {
      delay -> replaceSymbol (*par_name[q], *par_idx[q], &delay);
    }
    // replace time
    delay -> replaceSymbol (*time[0], *time_var, &delay);
    // take over ownership
    delayExpr[r].fromNode (delay);
  }

  // render delay into identity
  NodeEquals* del_rem = new NodeEquals;
  NodeFunction* del_id = new NodeFunction("delay", 2);
  del_id->addArgument (0, new NodeSymbol("a"));
  del_id->addArgument (1, new NodeSymbol("b"));
  del_rem->addArgument (0, del_id);
  del_rem->addArgument (1, new NodeSymbol("a"));

  //  replacing the variables
  for (size_t p = 0; p < var_name.size(); p++)
  {
    var_dot[p] -> replaceSymbol (*del_rem, &var_dot[p]);
    var_dot[p] -> replaceSymbol (*time[0], *time_var, &var_dot[p]);
    for (size_t q = 0; q < var_name.size(); q++)
    {
      var_dot[p] -> replaceSymbol (*var_name[q], *var_idx[q], &var_dot[p]);
    }
    for (size_t q = 0; q < par_name.size(); q++)
    {
      var_dot[p] -> replaceSymbol (*par_name[q], *par_idx[q], &var_dot[p]);
    }
  }
  del_rem -> deleteTree();
  delete del_rem;
  delete time_var;
  // END :: HANDLING THE DELAYS
    
  // fill up the RHS
  varDotExpr.resize (var_dot.size());
  for (size_t q = 0; q < var_dot.size(); q++)
  {
    varDotExpr[q].fromNode (var_dot[q]);
  }

  // cleaning up the variables
  for (auto it : var_name) { it->deleteTree(); delete it; }
  for (auto it : par_name) { it->deleteTree(); delete it; }
  for (auto it : par_value) { it->deleteTree(); delete it; }
  for (auto it : time) { it->deleteTree(); delete it; }
  for (auto it : macro) { it->deleteTree(); delete it; }
  for (auto it : var_idx) { it->deleteTree(); delete it; }
  for (auto it : par_idx) { it->deleteTree(); delete it; }
}

#include <mxml.h>

bool isValidName(const std::string s)
{
  const std::string lower("abcdefghijklmnopqrstuvwxyz");
  const std::string upper("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  const std::string letters = lower + upper + "_";
  const std::string digits("0123456789");
  const std::string validchars = letters + digits;

  if (s.empty())
  {
    return false;
  }
  if (letters.find(s[0]) == std::string::npos)
  {
    return false;
  }
  for (size_t k = 1; k < s.length(); ++k)
    if (validchars.find(s[k]) == std::string::npos)
    {
      return false;
    }
  return true;
}

void Expression::fromXML (std::string& oexpr, const std::string& xmlfilename)
{
  std::ostringstream out;
  FILE *xmlfile;
  mxml_node_t *tree;
  mxml_node_t *node;
  bool bad_attr;

  xmlfile = fopen(xmlfilename.c_str(), "r");
  P_ERROR_X3(xmlfile != NULL, "Error: Unable to open ", xmlfilename, "\n");
  tree = mxmlLoadFile(NULL, xmlfile, MXML_NO_CALLBACK);
  fclose(xmlfile);
  if (tree == NULL)
  {
    mxmlDelete(tree);
    P_MESSAGE3("Error: Unable to load the vector field from the file ", xmlfilename, 
               ". There may be an error in the XML definition of the vector field.");
  }

  node = mxmlFindElement(tree, tree, "VectorField", NULL, NULL, MXML_DESCEND);
  if (node == NULL)
  {
    mxmlDelete(tree);
    P_MESSAGE1("Error: No VectorField element found in XML defintion.\n");
  }
  else
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "IndependentVariable" && attr != "Description")
      {
        mxmlDelete(tree);
        P_MESSAGE2("Error: The VectorField element has an unknown attribute: ", attr);
      }
    }
    const char *attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: The VectorField element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The VectorField Name \"", attr, "\" is not valid.");
      }
      out << "vfname()=" << attr << ";\n";
    }
    attr = mxmlElementGetAttr(node, "Description");
    if (attr != NULL)
    {
      out << "vfdescription()=" << attr << ";\n";
    }
    attr = mxmlElementGetAttr(node, "IndependentVariable");
    if (attr == NULL)
    {
      out << "time()=t;\n";
    } else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The VectorField IndependentVariable \"", attr, "\" is not valid.");
      }
      out << "time()=" << attr << ";\n";
    }
  }

  //
  // Get the constants
  //
  for (node = mxmlFindElement(tree, tree, "Constant", NULL, NULL, MXML_DESCEND);
       node != NULL;
       node = mxmlFindElement(node, tree, "Constant", NULL, NULL, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Value" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: A Constant element has an unknown attribute: ", attr, "Valid Constant attributes are: Name, Value, Description, Latex.");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: A Constant element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Constant Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << name << "=";
      attr = mxmlElementGetAttr(node, "Value");
      if (attr == NULL)
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Constant element with Name=\"", name, "\" has no Value attribute.");
      }
      else
      {
        out << attr << ";\n";
      }
    }
  }

  //
  // Get the parameters
  //
  for (node = mxmlFindElement(tree, tree, "Parameter", NULL, NULL, MXML_DESCEND);
       node != NULL;
       node = mxmlFindElement(node, tree, "Parameter", NULL, NULL, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "DefaultValue" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: A Parameter element has an unknown attribute: ", attr, "Valid Parameter attributes are: Name, DefaultValue, Description, Latex.");
        bad_attr = true;
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: A Parameter element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Parameter Name \"", attr, "\" is not valid.\n");
      }
      std::string name(attr);
      std::string descr;
//       Parameter *p = new Parameter(name);
      attr = mxmlElementGetAttr(node, "Description");
      if (attr != NULL)
      {
        descr = attr;
        std::transform(descr.begin(), descr.end(), descr.begin(), ::tolower);
      }
      if (descr == "period")
        out << "period(" << name << ")=";
      else
        out << "par(" << name << ")=";
      
      attr = mxmlElementGetAttr(node, "DefaultValue");
      if (attr != NULL)
      {
        out << attr << ";\n";
      }
    }
  }

  //
  // Get the auxiliary expressions
  //
  for (node = mxmlFindElement(tree, tree, "Expression", NULL, NULL, MXML_DESCEND);
       node != NULL;
       node = mxmlFindElement(node, tree, "Expression", NULL, NULL, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Formula" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: An Expression element has an unknown attribute: ", attr, "Valid Expression attributes are: Name, Formula, Description, Latex");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: An Expression element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Expression Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << name << "=";
      attr = mxmlElementGetAttr(node, "Formula");
      if (attr == NULL)
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Expression with Name=\"", name, "\" has no Formula attribute.");
      }
      else
      {
        out << attr << ";\n";
      }
    }
  }

  //
  // Get the state variables
  //
  for (node = mxmlFindElement(tree, tree, "StateVariable", NULL, NULL, MXML_DESCEND);
       node != NULL;
       node = mxmlFindElement(node, tree, "StateVariable", NULL, NULL, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "DefaultInitialCondition" && attr != "Description"
          && attr != "Formula" && attr != "PeriodFrom" && attr != "PeriodTo"
          && attr != "DefaultHistory" && attr != "Mass" && attr != "Latex")
      {
      	mxmlDelete(tree);
        P_MESSAGE3("Error: A StateVariable element has an unknown attribute: ", attr, "Valid StateVariable attributes are: Name, Formula, Description, DefaultInitialCondition, DefaultHistory, PeriodFrom, PeriodTo, Latex.");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: A StateVariable element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The StateVariable Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << "dot(" << name << ")=";
      attr = mxmlElementGetAttr(node, "Formula");
      if (attr == NULL)
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The StateVariable with Name=\"", name, "\" has no Formula attribute.");
      }
      else
      {
        out << attr << ";\n";
      }
      attr = mxmlElementGetAttr(node, "DefaultInitialCondition");
      if (attr != NULL)
      {
        out << "init(" << name << ")=" << attr << ";\n";
      }
      attr = mxmlElementGetAttr(node, "Mass");
      if (attr != NULL)
      {
        out << "mass(" << name << ")=" << attr << ";\n";
      }
    }
  }

  //
  // Get the functions
  //
  for (node = mxmlFindElement(tree, tree, "Function", NULL, NULL, MXML_DESCEND);
       node != NULL;
       node = mxmlFindElement(node, tree, "Function", NULL, NULL, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Formula" && attr != "Description")
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: A Function element has an unknown attribute: ", attr, "Valid Function attributes are: Name, Formula, Description.");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == NULL)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: A Function element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Function Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << "fun(" << name << ")=";
      attr = mxmlElementGetAttr(node, "Formula");
      if (attr == NULL)
      {
        mxmlDelete(tree);
        P_MESSAGE3("Error: The Function element with Name=\"", name, "\" has no Formula attibute.\n");
      }
      else
      {
        out << attr << ";\n";
      }
    }
  }

  mxmlDelete(tree);
  oexpr = out.str ();
}

} // namespace
