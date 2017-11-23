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

#ifdef DEBUG
std::list<const Node *> Node::instances;
#endif

void Node::find (TokenType tp, const std::function<void(const Node*)>& found) const
{
  if (tp == type) found (this);
  for (size_t k = 0; k < children(); k++) getChild (k) -> find (tp, found);
}

const std::function<bool(const Node*)> Node::alwaysFalse = [] (const Node*) { return false; };

struct TokenPrecedence
{
  int precedence;
  bool leftAssoc;
  size_t nargs;
};

class Token
{
  public:
    Token (TokenType t, const std::string& n, const size_t loc) : type(t), name(n), vfloc(loc) { }
    bool isId() const { return type == TokenType::Symbol; }
    bool isOp() const { return (type != TokenType::Invalid)&&(type != TokenType::LeftParen)&&(type != TokenType::RightParen)&&(type != TokenType::Symbol)&&(type != TokenType::Number)&&(type != TokenType::Comma); }
    bool isSym() const { return type == TokenType::Symbol; }
    bool isLeft() const { return type == TokenType::LeftParen; }
    bool isRight() const { return type == TokenType::RightParen; }
// variables
    TokenType type;
    std::string name;
    const size_t vfloc;
};

class NodePlus : public Node
{
public:
  NodePlus (const size_t loc) : Node(TokenType::UnaryPlus, loc), arg(nullptr) {}
  void deleteTree () override { arg->deleteTree (); delete arg; }
  Node* copy () const override
  {
    NodePlus* cp = new NodePlus (vfloc);
    cp->arg = arg->copy();
    return cp;
  }
  size_t children () const override { return 1; }
  const Node* getChild (size_t n) const override { return arg; }
  Node* getChild (size_t n) override { return arg; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override
  { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const override
  {
    PE_ERROR_X1( arg != nullptr, vfloc, "Incerdibly wrong!");
    out << "(";
    out << "+";
    arg->print (out);
    out << ")";
  }
  bool operator != (const Node& node) const override
  {
    const NodePlus* nptr = dynamic_cast<const NodePlus*>(&node);
    if (nptr) return (*arg != *(nptr->arg));
    else return true;
  }
  bool operator < (const Node& node) const override
  {
    if (this->type < node.type) return true;
    else if (this->type == node.type)
    {
      if (*arg < *(static_cast<const NodePlus* >(&node)->arg)) return true;
      else return false;
    }
    return false;
  }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override
  {
    arg->replaceSymbol (sym, node, &arg);
  }
  void replaceSymbol (const Node& node, Node** parent) override
  {
    arg->replaceSymbol (node, &arg);
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override
  {
    arg->findFunction (lst, nm);
  }
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override { return arg->derivative (var, zero); }
  Node* optimize (const std::function<bool(const Node*)>& zero) override
  {
    arg = node_optimize (arg, zero);
    return arg;
  }
  // specific
  void addArgument (size_t pos, Node* arg_) { arg = arg_; }
private:
  Node* arg;
};

class NodeEquals : public Node
{
private:
  NodeEquals (TokenType tp, const size_t loc) : Node(tp, loc), args(2) {}
public:
  NodeEquals (const size_t loc) : Node(TokenType::Equals, loc), args(2) {}
  NodeEquals (bool, const size_t loc) : Node(TokenType::Define, loc), args(2) {}
  void deleteTree () override { for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; } }
  Node* copy () const override
  {
    NodeEquals* cp = new NodeEquals (type, vfloc);
    for(size_t k=0; k<args.size(); k++)
    {
      cp->args[k] = args[k]->copy();
    }
    return cp;
  }
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override { return args[n]; }
  Node* getChild (size_t n) override { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override
  { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const override
  {
    PE_ERROR_X1( args[0] != nullptr, vfloc, "Incerdibly wrong!");
    PE_ERROR_X1( args[1] != nullptr, vfloc, "Incerdibly wrong!");
    args[0]->print (out);
    if (type == TokenType::Equals) out << "=";
    else out << ":=";
    args[1]->print (out);
  }
  bool operator != (const Node& node) const override
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
  bool operator < (const Node& node) const override
  { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED!\n"); return true; }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override
  {
    for (size_t k=0; k<args.size(); k++) args[k]->replaceSymbol (sym, node, &args[k]);
  }
  void replaceSymbol (const Node& node, Node** parent) override
  {
    for (size_t k=0; k<args.size(); k++)
    {
      PE_ERROR_X1 (args[k] != nullptr, vfloc, "terribly wrong");
      args[k]->replaceSymbol (node, &args[k]);
    }
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override
  {
    for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
  }
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override { return new NodeNumber (0.0, vfloc); }
  Node* optimize (const std::function<bool(const Node*)>& zero) override
  {
    for(size_t k=0; k<args.size(); k++)
    {
      args[k] = node_optimize (args[k], zero);
    }
    return this;
  }
  // specific
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
private:
  std::vector<Node*> args;
};

class NodeSemicolon : public Node
{
public:
  NodeSemicolon (const size_t loc) : Node (TokenType::Semicolon, loc) {}
  void deleteTree () override { for(auto k=args.begin(); k!=args.end(); k++) { (*k)->deleteTree (); delete *k; } }
  Node* copy () const override
  {
    NodeSemicolon* cp = new NodeSemicolon (vfloc);
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      cp->args.push_back((*k)->copy());
    }
    return cp;
  }
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override
  {
    size_t j=0;
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      if (j == n) return *k;
      j++;
    }
    return nullptr;
  }
  Node* getChild (size_t n) override
  {
    size_t j=0;
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      if (j == n) return *k;
      j++;
    }
    return nullptr;
  }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override
  { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED\n"); return 0; }
  void print (std::ostream& out) const override
  {
    for(auto k=args.begin(); k!=args.end(); k++)
    {
      PE_ERROR_X1( (*k) != nullptr, vfloc, "Incerdibly wrong!");
      (*k)->print (out);
      out << ";\n";
    }
  }
  bool operator != (const Node& node) const override { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED!\n"); return true; }
  bool operator < (const Node& node) const override { PE_MESSAGE1(vfloc, "NOT IMPLEMENTED!\n"); return true; }
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->replaceSymbol (sym, node, &(*k));
  }
  void replaceSymbol (const Node& node, Node** parent) override
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->replaceSymbol (node, &(*k));
  }
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override
  {
    for(auto k=args.begin(); k!=args.end(); k++) (*k)->findFunction (lst, nm);
  }
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override { return new NodeNumber (0.0, vfloc); }
  // specific
  Node* optimize (const std::function<bool(const Node*)>& zero) override;
  void addArgumentBack (Node* arg) { args.push_back (arg); }
  void addArgumentFront (Node* arg) { args.push_front (arg); }
  void toList(std::list<Node*>& lst)
  {
    lst = std::move(args);
  }
private:
  std::list<Node*> args;
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
  const NodeSymbol* symref = dynamic_cast<const NodeSymbol*>(&sym);
  if (symref != nullptr)
  {
    if (name == symref->name)
    {
      Node* cp = node.copy();
      (*parent)->deleteTree ();
      P_ASSERT_X1 ((*parent) == this, "Wrong parent node.");
      delete *parent;
      *parent = cp;
    }
  }
}

void NodeSymbol::replaceSymbol (const Node& node, Node** parent)
{
  const NodeEquals* nd = dynamic_cast<const NodeEquals*>(&node);
  if (nd != nullptr)
  {
    const NodeSymbol* left = dynamic_cast<const NodeSymbol*>(nd->getChild(0));
    const Node* right = nd->getChild(1);
    if (left != nullptr)
    {
      replaceSymbol (*left, *right, parent);
    }
  } else PE_MESSAGE1(vfloc, "Not a definition\n");
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
  const NodeVar* symref = dynamic_cast<const NodeVar*>(&sym);
  if (symref)
  {
    if (idx == symref->idx)
    {
      Node* cp = node.copy ();
      (*parent)->deleteTree ();
      P_ASSERT_X1 ((*parent) == this, "Wrong parent node.");
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
  const NodePar* symref = dynamic_cast<const NodePar*>(&sym);
  if (symref)
  {
    if (idx == symref->idx)
    {
      Node* cp = node.copy ();
      (*parent)->deleteTree ();
      P_ASSERT_X1 ((*parent) == this, "Wrong parent node.");
      delete *parent;
      *parent = cp;
    }
  }
}

// NodeExpr
bool NodeExpr::operator != (const Node& node) const
{
  const NodeExpr* nptr = dynamic_cast<const NodeExpr*>(&node);
  if (nptr ==  nullptr) return true;
  else if (idx != nptr->idx) return true;
  else if (np != nptr->np) return true;
  else if (nv != nptr->nv) return true;
  else if ((np == 1) && (dp != nptr->dp)) return true;
  else if ((nv > 0) && (dv[0] != nptr->dv[0])) return true;
  else if ((nv > 1) && (dv[1] != nptr->dv[1])) return true;
  else return false;
}

bool NodeExpr::operator < (const Node& node) const
{
  if (type < node.type) return true;
  else if (type == node.type)
  {
    const NodeExpr* ne = static_cast<const NodeExpr* >(&node);
    if (idx < ne->idx) return true;
    else if (idx == ne->idx)
    {
      if (np < ne->np) return true;
      else if (np == ne->np)
      {
        if (nv < ne->nv) return true;
        else if (nv == ne->nv)
        {
          if ((np == 1) && (dp < ne->dp)) return true;
          else if ((np == 1) && (dp == ne->dp))
          {
            if ((nv > 0) && (dv[0] < ne->dv[0])) return true;
            else if ((nv > 0) && (dv[0] == ne->dv[0]))
            {
              if ((nv > 1) && (dv[1] < ne->dv[1])) return true;
              else return false;
            }
          }
        }
      }
    }
  }
  return false;
}

void NodeExpr::replaceSymbol (const Node& sym, const Node& node, Node** parent)
{
  if (this->operator==(sym))
  {
    Node* cp = node.copy ();
    (*parent)->deleteTree ();
    P_ASSERT_X1 ((*parent) == this, "Wrong parent node.");
    delete *parent;
    *parent = cp;
  }
}

// NodeFunction
Node* NodeFunction::copy () const
{
  NodeFunction* cp = new NodeFunction (name, args.size(), vfloc);
  for(size_t k=0; k<args.size(); k++) cp->args[k] = args[k]->copy();
  return cp;
}

void NodeFunction::print (std::ostream& out) const
{
  if (name == "abs") out << "std::fabs(";
  else if (name == "sign") out << "sign(";
  else if (name == "heaviside") out << "heaviside(";
  else out << "std::" << name << "(";
  for(size_t k=0; k<args.size(); k++)
  {
    PE_ERROR_X1( args[k] != nullptr, vfloc, "Incerdibly wrong!");
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
    const NodeFunction* macro = dynamic_cast<const NodeFunction*>(node.getChild(0));
    const Node* right = node.getChild(1);
    if (macro != nullptr)
    {
      if ( (name == macro->name)&&(children() == macro->children()) )
      {
        Node* cp = right->copy ();
        for (size_t k=0; k<args.size(); k++)
        {
          if (macro->args[k]->type == TokenType::Symbol)
          {
            cp->replaceSymbol (*macro->args[k], *args[k], &cp);
          }
          else PE_MESSAGE1(vfloc, "Not a symbol");
        }
        (*parent)->deleteTree ();
        P_ASSERT_X1 ((*parent) == this, "Wrong parent node.");
        delete *parent;
        *parent = cp;
        return;
      }
    }
  }
}

void NodeFunction::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
  if (name == nm) lst.push_back (this);
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
    NodeFunctionA1 (const std::string& nm, const size_t loc);
    NodeFunctionA1 (const size_t loc);
    const NodeFunctionA1& operator = (const NodeFunctionA1&);
    Node* copy () const override;
    size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  private:
    std::function<double(double)> function;
};

NodeFunctionA1::NodeFunctionA1 (const std::string& nm, const size_t loc) : NodeFunction (nm, 1, loc)
{
  if (nm == "abs") function = [](double a) { return fabs(a); };
  else if (nm == "exp") function = [](double a) { return exp(a); };
  else if (nm == "log") function = [](double a) { return log(a); };
  else if (nm == "sin") function = [](double a) { return sin(a); };
  else if (nm == "cos") function = [](double a) { return cos(a); };
  else if (nm == "tan") function = [](double a) { return tan(a); };
  else if (nm == "asin") function = [](double a) { return asin(a); };
  else if (nm == "acos") function = [](double a) { return acos(a); };
  else if (nm == "atan") function = [](double a) { return atan(a); };
  else if (nm == "sinh") function = [](double a) { return sinh(a); };
  else if (nm == "cosh") function = [](double a) { return cosh(a); };
  else if (nm == "tanh") function = [](double a) { return tanh(a); };
  else if (nm == "asinh") function = [](double a) { return asinh(a); };
  else if (nm == "acosh") function = [](double a) { return acosh(a); };
  else if (nm == "atanh") function = [](double a) { return atanh(a); };
  else if (nm == "sign") function = &localSign;
  else if (nm == "heaviside") function = &localHeaviside;
  else function = nullptr;
}

Node* NodeFunctionA1::copy () const
{
  NodeFunctionA1* cp = new NodeFunctionA1 (NodeFunction::name, vfloc);
  for(size_t k=0; k<args.size(); k++) cp->args[k] = args[k]->copy();
  cp->function = function;
  return cp;
}

// NodeFunctionA2
class NodeFunctionA2 : public NodeFunction
{
  public:
    NodeFunctionA2 (const std::string& nm, const size_t loc);
    Node* copy () const override;
    size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  private:
    std::function<double(double,double)> function;
};

NodeFunctionA2::NodeFunctionA2 (const std::string& nm, const size_t loc) : NodeFunction (nm, 2, loc)
{
  if (nm == "pow") function = [](double a, double b) { return pow(a, b); };
}

Node* NodeFunctionA2::copy () const
{
  NodeFunctionA2* cp = new NodeFunctionA2 (NodeFunction::name, vfloc);
  for(size_t k=0; k<this->args.size(); k++) cp->args[k] = this->args[k]->copy();
  cp->function = function;
  return cp;
}

// NodeAdd

Node* NodeAdd::copy () const
{
  NodeAdd* cp = new NodeAdd (args.size(), vfloc);
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->copy();
    cp->mul[k] = mul[k];
  }
  return cp;
}

void NodeAdd::print (std::ostream& out) const
{
  PE_ERROR_X1 (args.size() > 0, vfloc, "terribly wrong");
  out << "(";
  for (size_t k=0; k<args.size(); k++)
  {
    PE_ERROR_X1( args[k] != nullptr, vfloc, "Incerdibly wrong!");
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
    PE_ERROR_X1 (args[k] != nullptr, vfloc, "terribly wrong");
    args[k]->replaceSymbol (node, &args[k]);
  }
}

void NodeAdd::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
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
  NodeTimes* cp = new NodeTimes (args.size(), vfloc);
  cp->smul = smul;
  for(size_t k=0; k<args.size(); k++)
  {
    PE_ERROR_X1 (args[k] != nullptr, vfloc, "terribly wrong");
    cp->args[k] = args[k] -> copy ();
    cp->divide[k] = divide[k];
  }
  return cp;
}

void NodeTimes::print (std::ostream& out) const
{
  PE_ERROR_X1 (args.size() > 0, vfloc, "terribly wrong");
  out << "(";
  if (smul != 1.0)
  {
    if (smul < 0.0) out << "-";
    if (smul != -1.0) out << fabs(smul) << "*";
  }
  for(size_t k=0; k<args.size(); k++)
  {
    PE_ERROR_X1 (args[k] != nullptr, vfloc, "terribly wrong");
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
    PE_ERROR_X4 (args[k] != nullptr, vfloc, "terribly wrong: ", k, " of ", args.size() );
    args[k]->replaceSymbol (node, &args[k]);
  }
}

void NodeTimes::findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const
{
  for (size_t k=0; k<args.size(); k++) args[k]->findFunction (lst, nm);
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
  NodePower* cp = new NodePower (vfloc);
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->copy();
  }
  return cp;
}

void NodePower::print (std::ostream& out) const
{
  out << "pow(";
  PE_ERROR_X1( args[0] != nullptr, vfloc, "Incerdibly wrong!");
  PE_ERROR_X1( args[1] != nullptr, vfloc, "Incerdibly wrong!");
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

Node* NodePower::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
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
  NodeTimes* t1 = new NodeTimes (2, vfloc);
  NodeTimes* t2 = new NodeTimes (3, vfloc);
  NodeTimes* t3 = new NodeTimes (2, vfloc);
  NodeAdd* s1 = new NodeAdd (2, vfloc);
  NodeAdd* s2 = new NodeAdd (2, vfloc);
  NodeFunctionA1* fun = new NodeFunctionA1 ("log", vfloc);
  NodePower* pw = new NodePower (vfloc);

  s1 -> addArgument (0, b -> copy (), 1.0);
  s1 -> addArgument (1, new NodeNumber (1.0, vfloc), -1.0);
  pw -> args[0] = a -> copy();
  pw -> args[1] = s1;
  t1 -> addArgument (0, a -> derivative (var, zero), 1.0, false);
  t1 -> addArgument (1, b -> copy (), 1.0, false);
  fun -> addArgument(0, a -> copy ());
  t2 -> addArgument (0, a -> copy (), 1.0, false);
  // Need to check if D[b] of t2 is zero -> log(negative) = NaN
  Node* tmpDb = b -> derivative (var, zero);
  Node* tmpLog = fun;

  tmpDb = node_optimize (tmpDb, Node::alwaysFalse);

  if (tmpDb->type == TokenType::Number)
  {
    const double v = static_cast<const NodeNumber*>(tmpDb)->getValue();
    if (v == 0)
    {
      fun->deleteTree ();
      delete fun;
      tmpLog = new NodeNumber(0.0, vfloc);
    }
  }
  t2 -> addArgument (1, tmpDb, 1.0, false);
  t2 -> addArgument (2, tmpLog, 1.0, false);
  s2 -> addArgument (0, t1, 1.0);
  s2 -> addArgument (1, t2, 1.0);
  t3 -> addArgument (0, pw, 1.0, false);
  t3 -> addArgument (1, s2, 1.0, false);
  Node* t4 = t3;
  t4 = node_optimize (t4, Node::alwaysFalse);
  return t4;
}

Node* NodePower::optimize (const std::function<bool(const Node*)>& zero)
{
  for (size_t k=0; k<args.size(); k++)
  {
    args[k] = node_optimize (args[k], zero);
  }
  // if both of them are numbers ...
  NodeNumber* a0 = dynamic_cast<NodeNumber*>(args[0]);
  NodeNumber* a1 = dynamic_cast<NodeNumber*>(args[1]);
  if ((a0 != nullptr)&&(a1 != nullptr))
  {
    NodeNumber* res = new NodeNumber (pow (a0->getValue(), a1->getValue()), vfloc);
    // should delete both arguments
    deleteTree ();
    return res;
  }
  if (a1 != nullptr)
  {
    if (a1 -> getValue () == 1.0)
    {
      Node* tmp = args[0];
      delete args[1];
      return tmp;
    }
    if (a1 -> getValue () == 0.0)
    {
      NodeNumber* res = new NodeNumber (1.0, vfloc);
      // should delete both arguments
      deleteTree ();
      return res;
    }
  }
  return this;
}

NodeFunction* NodeFunction::specialize ()
{
  if (args.size() == 1)
  {
    NodeFunctionA1* fun = new NodeFunctionA1 (name, vfloc);
    fun -> addArgument (0, args[0]);
    return fun;
  }
  if (args.size() == 2)
  {
    NodeFunctionA2* fun = new NodeFunctionA2 (name, vfloc);
    fun -> addArgument (0, args[0]);
    fun -> addArgument (1, args[1]);
    return fun;
  }
  return this;
}

Node* NodeFunction::optimize (const std::function<bool(const Node*)>& zero)
{
  for (size_t k=0; k<args.size(); k++)
  {
    args[k] = node_optimize (args[k], zero);
  }
  // replace power with NodePower
  if (name == "pow")
  {
    NodePower* pw = new NodePower (vfloc);
    pw -> addArgument (0, args[0]);
    pw -> addArgument (1, args[1]);
    return pw;
  }
  bool isnum = true;
  for (size_t k=0; k<args.size(); k++)
  {
    if (args[k]->type != TokenType::Number) isnum = false;
  }
  // if all are numbers -> evaluate
  if (isnum)
  {
    std::vector<double> par;
    double res;
    this->evaluateWithPar (&res, par);
    NodeNumber* num = new NodeNumber (res, vfloc);
    return num;
  }
  // this is true all the time!
  return specialize ();
}

Node* NodeAdd::optimize (const std::function<bool(const Node*)>& zero)
{
  double sum = 0;   // the constant part
  std::vector<Node*> newargs;
  std::vector<double> newmul;
  for (size_t k = 0; k < args.size(); k++)
  {
    args[k] = node_optimize (args[k], zero);
  }
  for (size_t k = 0; k < args.size(); k++)
  {
    if (mul[k] == 0.0)
    {
      args[k]->deleteTree ();
      delete args[k];
      args[k] = nullptr;
      continue;
    }
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[k]);
    if (num)
    {
      sum += mul[k] * (num -> getValue());
      args[k]->deleteTree ();
      delete args[k];
      args[k] = nullptr;
      continue;
    }
    NodeAdd* child = dynamic_cast<NodeAdd*>(args[k]);
    if (child)
    {
      for (size_t p = 0; p < child->args.size(); p++)
      {
        newargs.push_back (child->args[p]);
        newmul.push_back (mul[k]*child->mul[p]);
      }
      delete child;
      args[k] = nullptr;
      continue;
    }
    NodeTimes* times = dynamic_cast<NodeTimes*>(args[k]);
    if (times)
    {
      newargs.push_back (args[k]);
      newmul.push_back (mul[k]*(times->getMul()));
      times->setMul(1.0);
      args[k] = nullptr;
      continue;
    }
    newargs.push_back (args[k]);
    args[k] = nullptr;
    newmul.push_back (mul[k]);
  }
  args = newargs;
  mul = newmul;
  if (sum != 0.0)
  {
    args.push_back (new NodeNumber (sum, vfloc));
    mul.push_back (1.0);
  }
  if (args.size() == 0)
  {
    Node* tmp = new NodeNumber (0.0, vfloc);
    return tmp;
  }
  if (args.size() == 1)
  {
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[0]);
    if (num)
    {
      Node* tmp = new NodeNumber (mul[0] * (num ->getValue()), vfloc);
      deleteTree();
      return tmp;
    } else
    {
      if (mul[0] == 1.0)
      {
        // taking over ownership
        Node* tmp = args[0];
        args[0] = nullptr;
        return tmp;
      } else
      {
        // convert into times
        NodeTimes* tmp = new NodeTimes (1, vfloc);
        tmp -> addArgument (0, args[0], mul[0], false);
        args[0] = nullptr;
        return tmp;
      }
    }
  }
  return this;
}

Node* NodeTimes::optimize (const std::function<bool(const Node*)>& zero)
{
  std::vector<Node*> newargs;
  std::vector<bool> newdivide;
  double newsmul = smul;
  for (size_t k = 0; k < args.size(); k++)
  {
    P_ASSERT_X1 (args[k] != nullptr, "Optimizing null pointer.");
    args[k] = node_optimize (args[k], zero);
  }
  for (size_t k = 0; k < args.size(); k++)
  {
//     args[k]->print (std::cout); std::cout << "\t";
    NodeNumber* num = dynamic_cast<NodeNumber*>(args[k]);
    if (num)
    {
      if (divide[k]) newsmul /= num -> getValue();
      else newsmul *= num -> getValue();
      // deleting the the term
      args[k]->deleteTree ();
      delete args[k];
      args[k] = nullptr;
      continue;
    }
    NodeTimes* child = dynamic_cast<NodeTimes*>(args[k]);
    if (child)
    {
      if (divide[k]) newsmul /= child->smul;
      else newsmul *= child->smul;
      for (size_t p = 0; p < child->children(); p++)
      {
        // taking over ownership
        newargs.push_back (child->args[p]);
        child->args[p] = nullptr;
        if (divide[k]) newdivide.push_back (!child->divide[p]);
        else newdivide.push_back (child->divide[p]);
      }
      delete args[k];
      args[k] = nullptr;
      continue;
    }
    if (args[k])
    {
      newargs.push_back (args[k]);
      args[k] = nullptr;
      newdivide.push_back (divide[k]);
    }
  }
  args = newargs;
  divide = newdivide;
  smul = newsmul;
  if (smul == 0.0)
  {
    deleteTree ();
    return new NodeNumber (0.0, vfloc);
  }
  // single argument
  if ((smul == 1.0)&&(args.size() == 1)&&(!divide[0]))
  {
    // self destruct
    Node* tmp = args[0];
    return tmp;
  }
  if (args.size() == 0)
  {
    const double smul_copy = smul;
    return new NodeNumber (smul_copy, vfloc);
  }
  return this;
}

Node* NodeSemicolon::optimize (const std::function<bool(const Node*)>& zero)
{
  auto it=args.begin();
  while (it != args.end())
  {
    (*it) = node_optimize (*it, zero);
    NodeSemicolon* child = dynamic_cast<NodeSemicolon*>(*it);
    if (child)
    {
      for (auto p=child->args.begin(); p != child->args.end(); p++)
      {
        (*p) = node_optimize (*p, zero);
        args.insert (it, *p);
      }
      it = args.erase (it);
      delete child;
    } else
    it++;
  }
  return this;
}

// DERIVATIVE IMPLEMENTATIONS

Node* NodeNumber::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  return new NodeNumber (0.0, vfloc);
}

Node* NodeSymbol::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  if (var->type == TokenType::Symbol)
  {
    const NodeSymbol* node = static_cast<const NodeSymbol* >(var);
    {
      if (node->name == name)
      {
        return new NodeNumber (1.0, vfloc);
      }
    }
  }
  return new NodeNumber (0.0, vfloc);
}

Node* NodeVar::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  if (var->type == TokenType::Var)
  {
    const NodeVar* node = static_cast<const NodeVar* >(var);
    {
      if (node->idx == idx)
      {
        return new NodeNumber (1.0, vfloc);
      }
    }
  }
  return new NodeNumber (0.0, vfloc);
}

Node* NodePar::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  if (var->type == TokenType::Par)
  {
    const NodePar* node = static_cast<const NodePar* >(var);
    {
      if (node->idx == idx)
      {
        return new NodeNumber (1.0, vfloc);
      }
    }
  }
  return new NodeNumber (0.0, vfloc);
}

Node* NodeExpr::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  PE_ERROR_X1 (nv + np < 2, vfloc, "Up to second order derivatives are allowed.");
  const NodeVar* nvp = dynamic_cast<const NodeVar* >(var);
  if (nvp != nullptr)
  {
    NodeExpr* ne = this->copy();
    ne->dv[nv] = nvp->getIdx();
    ne->nv += 1;
    if (zero (ne))
    {
      delete ne;
      return new NodeNumber (0.0, vfloc);
    } else
    {
      return ne;
    }
  }
  const NodePar* npp = dynamic_cast<const NodePar* >(var);
  if (npp != nullptr)
  {
    PE_ERROR_X1 (np < 1, vfloc, "Was already differentiated w.r.t. a parameter.");
    NodeExpr* ne = this->copy();
    ne->dp = npp->getIdx();
    ne->np += 1;
    if (zero (ne))
    {
      delete ne;
      return new NodeNumber (0.0, vfloc);
    } else
    {
      return ne;
    }
  }
  return new NodeNumber (0.0, vfloc);
}

Node* NodeAdd::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  NodeAdd* cp = new NodeAdd (args.size(), vfloc);
  for(size_t k=0; k<args.size(); k++)
  {
    cp->args[k] = args[k]->derivative (var, zero);
    cp->mul[k] = mul[k];
  }
  return cp;
}

Node* NodeTimes::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  // f = (f1 f2 f2 ..) / (g1 g2 g3 ..)
  NodeAdd* root = new NodeAdd (args.size(), vfloc);
  for(size_t k=0; k<args.size(); k++)
  {
    NodeTimes* cp = new NodeTimes (0, vfloc);
    cp->smul = smul;
    if (divide[k])
    {
      // divisions
      cp->smul *= -1.0;
      cp->args.push_back (args[k]->derivative (var, zero));
      cp->divide.push_back (false);
      cp->args.push_back (args[k]->copy ());
      cp->divide.push_back (true);
      cp->args.push_back (args[k]->copy ());
      cp->divide.push_back (true);
    }
    else
    {
      P_ASSERT_X1(args[k] != nullptr, "Nullptr\n");
      Node* d = args[k]->derivative (var, zero);
      cp->args.push_back (d);
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
  Node* res = root;
  res = node_optimize (res, Node::alwaysFalse);
  return res;
}

Node* NodeFunction::derivative (const Node* var, const std::function<bool(const Node*)>& zero) const
{
  if (name == "abs")
  {
    // -> sign
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("sign", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "exp")
  {
    // -> exp
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("exp", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "log")
  {
    // -> 1/x
    NodeTimes* root = new NodeTimes (2, vfloc);
    root -> addArgument (0, args[0] -> copy (), 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "sin")
  {
    // -> cos
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("cos", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "cos")
  {
    // -> -sin
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("sin", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, -1.0, false);
  root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
  return root;
  }
  else if (name == "tan")
  {
    // -> cos^(-2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("cos", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    NodePower* pow = new NodePower (vfloc);
    pow -> addArgument (0, fun);
    pow -> addArgument (1, new NodeNumber (2.0, vfloc));
    root -> addArgument (0, pow, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "asin")
  {
    // -> 1/{(1-x^2)^(1/2)}
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (0, new NodeNumber (1.0, vfloc), 1.0);
    add -> addArgument (1, pow_square, -1.0);
    NodePower* pow_sqrt = new NodePower (vfloc);
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber (0.5, vfloc));
    root -> addArgument (0, pow_sqrt, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "acos")
  {
    // -> -1/{(1-x^2)^(1/2)}
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (0, new NodeNumber (1.0, vfloc), 1.0);
    add -> addArgument (1, pow_square, -1.0);
    NodePower* pow_sqrt = new NodePower (vfloc);
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber (0.5, vfloc));
    root -> addArgument (0, pow_sqrt, -1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "atan")
  {
    // -> 1/(1-x^2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (0, new NodeNumber (1.0, vfloc), 1.0);
    add -> addArgument (1, pow_square, -1.0);
    root -> addArgument (0, add, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "sinh")
  {
    // -> cosh
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("cosh", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "cosh")
  {
    // -> sinh
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("sinh", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    root -> addArgument (0, fun, 1.0, false);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "tanh")
  {
    // -> cosh^(-2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow = new NodePower (vfloc);
    NodeFunctionA1* fun = new NodeFunctionA1 ("cosh", vfloc);
    fun -> addArgument (0, args[0] -> copy ());
    pow -> addArgument (0, fun);
    pow -> addArgument (1, new NodeNumber (2.0, vfloc));
    root -> addArgument (0, pow, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "asinh")
  {
    // -> 1/(1+x^2)^(1/2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (0, new NodeNumber (1.0, vfloc), 1.0);
    add -> addArgument (1, pow_square, 1.0);
    NodePower* pow_sqrt = new NodePower (vfloc);
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber (0.5, vfloc));
    root -> addArgument (0, pow_sqrt, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "acosh")
  {
    // -> 1/(x^2-1)^(1/2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (1, new NodeNumber (1.0, vfloc), -1.0);
    add -> addArgument (0, pow_square, 1.0);
    NodePower* pow_sqrt = new NodePower (vfloc);
    pow_sqrt -> addArgument (0, add);
    pow_sqrt -> addArgument (1, new NodeNumber (0.5, vfloc));
    root -> addArgument (0, pow_sqrt, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "atanh")
  {
    // -> 1/(1-x^2)
    NodeTimes* root = new NodeTimes (2, vfloc);
    NodePower* pow_square = new NodePower (vfloc);
    pow_square -> addArgument (0, args[0] -> copy ());
    pow_square -> addArgument (1, new NodeNumber (2.0, vfloc));
    NodeAdd* add = new NodeAdd (2, vfloc);
    add -> addArgument (0, new NodeNumber (1.0, vfloc), 1.0);
    add -> addArgument (1, pow_square, -1.0);
    root -> addArgument (0, add, 1.0, true);
    root -> addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
    return root;
  }
  else if (name == "pow")
  {
    // promote to NodePower
    NodePower* pw = new NodePower (vfloc);
    pw -> addArgument (0, args[0]);
    pw -> addArgument (1, args[1]);
    Node* retval = pw -> derivative (var, zero);
    delete pw;
    return retval;
  }
  else if (name == "sign")
  {
    return new NodeNumber (0.0, vfloc);
  }
  else if (name == "heaviside")
  {
    return new NodeNumber (0.0, vfloc);
  }
  else
  {
    // this is when we don't know the function but do not want to leave the root hanging
    // -> cos
    NodeAdd* add = new NodeAdd (args.size(), vfloc);
    for (size_t k=0; k<args.size(); ++k)
    {
      std::ostringstream str;
      str << name << "_d" << k;
      NodeFunction* fun = new NodeFunction (str.str (), args.size(), vfloc);
      NodeTimes* times = new NodeTimes (2, vfloc);
      add->addArgument(k, times, 1.0);
      times->addArgument(0, fun, 1.0, false);
      times->addArgument (1, args[0] -> derivative (var, zero), 1.0, false);
      for (size_t p=0; p<args.size(); ++p)
      {
        fun -> addArgument (p, args[p] -> copy ());
      }
    }
    return add;
  }
}

// Evaluating the expression

void Node::evaluateWithPar (double *res, const std::vector<double>& parInit)
{
  ValueStack stk (res, 1);
  auto fun = [&stk, parInit] (size_t sp, const Node* node)
  {
    const NodePar* par_node = dynamic_cast<const NodePar*>(node);
    if (par_node != nullptr)
    {
      stk[sp].data[0] = parInit [par_node->getIdx()];
      stk[sp].skip = 1;
    }
  };
  evaluate (stk, 0, fun, 1);
}

// Vector AND Function arguments

size_t NodeNumber::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack should be pre-allocated
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = value;
  }
  return sp;
}

size_t NodeSymbol::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // this is an unresolved symbol, should not be called
  PE_MESSAGE2(vfloc, "Undefined symbol : ", name);
  // stack should be pre-allocated
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = 0;
  }
  return sp;
}

size_t NodeVar::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack should be pre-allocated
  fun (sp, this);
  return sp;
}

size_t NodePar::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack should be pre-allocated
  fun (sp, this);
  return sp;
}

size_t NodeExpr::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack should be pre-allocated
  fun (sp, this);
  return sp;
}

size_t NodeFunction::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  PE_MESSAGE3(vfloc, "NodeFunction::evaluate : it is a base class, cannot be called ", name, "\n");
  return sp;
}

size_t NodeFunctionA1::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack should be pre-allocated
  args[0] -> evaluate(stack, sp, fun, len);
  double* data = stack[sp].data;
  size_t skip = stack[sp].skip;
  if (!(this->function)) { PE_MESSAGE3 (vfloc, "'", name, "' is not defined as a function."); }
  for (size_t k = 0; k < len; k++)
  {
    data[skip*k] = this->function(data[skip*k]);
  }
  return sp;
}

size_t NodeFunctionA2::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack needs one more element
  if (stack.size() - 1 < sp + 1) stack.resizeDepth (sp + 2 + KNUT_STACK_INCREMENT);
  args[0] -> evaluate(stack, sp, fun, len);
  args[1] -> evaluate(stack, sp+1, fun, len);
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

size_t NodeAdd::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack needs more elements
  if (stack.size() - 1 < sp + args.size())
  {
    stack.resizeDepth (sp + args.size() + KNUT_STACK_INCREMENT);
  }
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, fun, len); }
  if (spi != sp + args.size()) PE_MESSAGE7(vfloc, "NodeAdd::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
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

size_t NodeTimes::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack needs more elements
  if (stack.size() - 1 < sp + args.size())
  {
    stack.resizeDepth (sp + args.size() + KNUT_STACK_INCREMENT);
  }
  size_t spi = sp;
  for (size_t k=0; k<args.size(); k++) { spi = args[k] -> evaluate(stack, sp+k+1, fun, len); }
  if (spi != sp + args.size()) PE_MESSAGE7(vfloc, "NodeTimes::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size(), "\n");
  for (size_t p = 0; p < len; p++)
  {
    stack[sp].data[p*stack[sp].skip] = smul;
  }
  for (size_t k = 0; k < args.size(); k++)
  {
    for (size_t p = 0; p < len; p++)
    {
      if (divide[k])
      {
        if ( stack[sp+k+1].data[p*stack[sp+k+1].skip] == 0 ) 
          PE_MESSAGE7(vfloc, "NodeTimes::evaluate : Division by zero ", spi, " ? ", sp, " + ", args.size(), "\n");
        stack[sp].data[p*stack[sp].skip] /= stack[sp+k+1].data[p*stack[sp+k+1].skip];
      }
      else stack[sp].data[p*stack[sp].skip] *= stack[sp+k+1].data[p*stack[sp+k+1].skip];
    }
  }
  return sp;
}

size_t NodePower::evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const
{
  // stack needs more elements
  if (stack.size() - 1 < sp + args.size())
  {
    stack.resizeDepth (sp + args.size() + KNUT_STACK_INCREMENT);
  }
  args[0] -> evaluate(stack, sp, fun, len);
  size_t spi = args[1] -> evaluate(stack, sp+1, fun, len);
  if (spi != sp + args.size() - 1) PE_MESSAGE7(vfloc, "NodePower::evaluate : Stack error ", spi, " ? ", sp, " + ", args.size()-1, "\n");
  for (size_t p = 0; p < len; p++)
  {
    stack[sp].data[p*stack[sp].skip] = pow (stack[sp].data[p*stack[sp].skip], stack[sp+1].data[p*stack[sp+1].skip]);
  }
  return sp;
}

// There are a number of limits to tokens
// alphanumeric is an identifier
// operators inbetween are the delimiters and they are tokens separately
//

enum class CharType : int { Op, Identifier, Space, Invalid };

static inline CharType findCharType (char c)
{
  if ( isalnum ( c ) || (c == '_') || (c == '.') ) return CharType::Identifier;
  else if ( (c == '(') || (c == ')') || (c == '*') || (c == '/') || (c == '+') || (c == '-') || (c == '^') || (c == ',') || (c == '=') || (c == ':') || (c == ';')) return CharType::Op;
  else if ( isspace ( c ) ) return CharType::Space;
  else return CharType::Invalid;
}

// This decomposes the program file into tokens
// TODO we need better error handling here
// need to seperate function recognition from identifier
void tokenize (std::vector<Token>& stream, const std::map<TokenType, TokenPrecedence>& token_table, const std::string& expression)
{
  const char *str = expression.c_str ();
  const size_t len = expression.size ();

  size_t k = 0;
  while (k < len)
  {
    // ignore comments
    if (str[k] == '#')
    {
      while ((str[k] != '\n')&&(str[k] != '\r')&&(str[k] != '\v')&&(str[k] != '\f')) k++;
      continue;
    }
    // ignore space
    if (findCharType (str[k]) == CharType::Space)
    {
      k++;
      continue;
    } else if (findCharType (str[k]) == CharType::Identifier)
    {
      if (std::isdigit(str[k]) || str[k] == '.')
      {
        char* endptr = nullptr;
        strtod (&str[k], &endptr);
        PE_ERROR_X1(endptr != &str[k], k, "Cannot convert number.");
        size_t p = (endptr - &str[k])/sizeof(char);
        std::string tn (&str[k], p);
        // TODO add to a symbol table and only use symbol ID later on
        stream.push_back (Token (TokenType::Symbol, tn, k) );
        k = k + p;
        continue;
      } else
      {
        size_t p = 0;
        while (findCharType (str[k+p]) == CharType::Identifier) p++;
        std::string tn (&str[k], p);
        // TODO add to a symbol table and only use symbol ID later on
        stream.push_back (Token (TokenType::Symbol, tn, k) );
        k = k + p;
        continue;
      }
    } else if (findCharType (str[k]) == CharType::Op)
    {
      std::string tn (&str[k],1);
      // a single character
      switch (tn[0])
      {
        case '+':
          stream.push_back (Token (TokenType::Plus, tn, k) );
          break;
        case '-':
          stream.push_back (Token (TokenType::Minus, tn, k) );
          break;
        case '*':
          stream.push_back (Token (TokenType::Times, tn, k) );
          break;
        case '/':
          stream.push_back (Token (TokenType::Divide, tn, k) );
          break;
        case '^':
          stream.push_back (Token (TokenType::Power, tn, k) );
          break;
        case '(':
          stream.push_back (Token (TokenType::LeftParen, tn, k) );
          break;
        case ')':
          stream.push_back (Token (TokenType::RightParen, tn, k) );
          break;
        case ',':
          stream.push_back (Token (TokenType::Comma, tn, k) );
          break;
        case '=':
          stream.push_back (Token (TokenType::Equals, tn, k) );
          break;
        case ':':
          k++;
          if (str[k] == '=') stream.push_back (Token (TokenType::Define, ":=", k-1) );
          else
          {
            k--;
            PE_MESSAGE3(k, "Unknown operator '", tn, "'.\n");
          }
          break;
        case ';':
          stream.push_back (Token (TokenType::Semicolon, tn, k) );
          break;
        default:
          PE_MESSAGE3(k, "Unknown operator '", tn, "'.\n");
      }
      k++;
    } else
    {
      PE_MESSAGE3(k, "Invalid character '", std::string(1,str[k]), "'\n");
      k++; // it is never reached, so this is dead code
    }
  }
  // this is now detecting token types : this should be moved
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
  else PE_MESSAGE1(stream.back().vfloc, "Missing semicolon\n");
}

// this converts tokens into nodes
void opToTree (std::list<Node*>& ts, const Token& tk, size_t nargs)
{
  if (ts.size() >= nargs)
  {
    switch (tk.type)
    {
      case TokenType::LeftParen:
      case TokenType::RightParen:
        PE_MESSAGE3(tk.vfloc, "Unexpected Parenthesis ", tk.name, ".");
        break;
      case TokenType::Symbol:
      {
        char* endptr;
        double res = strtod (tk.name.c_str(), &endptr);
        if (endptr != tk.name.c_str())
        {
          NodeNumber* node = new NodeNumber (res, tk.vfloc);
          ts.push_back (node);
        } else
        {
          NodeSymbol* node = new NodeSymbol (tk.name, tk.vfloc);
          ts.push_back (node);
        }
      }
      break;
      case TokenType::Function:
      {
        if (nargs == 1)
        {
          NodeFunctionA1* node = new NodeFunctionA1 (tk.name, tk.vfloc);
          // put the argument from the stack
          node->addArgument (0, ts.back());
          ts.pop_back ();
          ts.push_back (node);
        } else if (nargs == 2)
        {
          NodeFunctionA2* node = new NodeFunctionA2 (tk.name, tk.vfloc);
          // put the argument from the stack
          node->addArgument (1, ts.back());
          ts.pop_back ();
          node->addArgument (0, ts.back());
          ts.pop_back ();
          ts.push_back (node);
        } else
        {
          NodeFunction* node = new NodeFunction (tk.name, nargs, tk.vfloc);
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
        // Basically use ` + <value>'
        if (nargs == 1)
        {
          NodeAdd* node = new NodeAdd (1, tk.vfloc);
          node->addArgument (0, ts.back(), 1.0);
          ts.pop_back ();
          node->sort ();
          ts.push_back (node);
        } else PE_MESSAGE1(tk.vfloc, "Wrong number of arguments.");
      }
      break;
      case TokenType::UnaryMinus:
      {
        // Basically use ` - <value>'
        if (nargs == 1)
        {
          NodeAdd* node = new NodeAdd (1, tk.vfloc);
          node->addArgument (0, ts.back(), -1.0);
          ts.pop_back ();
          node->sort ();
          ts.push_back (node);
        } else PE_MESSAGE1(tk.vfloc, "Wrong number of arguments.");
      }
      break;
      case TokenType::Plus:
      {
        NodeAdd* node = new NodeAdd (nargs, tk.vfloc);
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
        NodeAdd* node = new NodeAdd (nargs, tk.vfloc);
        // put the argument from the stack
        if (nargs == 2)
        {
          node->addArgument (1, ts.back(), -1.0);
          ts.pop_back ();
          node->addArgument (0, ts.back(), 1.0);
          ts.pop_back ();
        } else PE_MESSAGE1(tk.vfloc, "Wrong number of arguments.");
        node->sort ();
        ts.push_back (node);
      }
      break;
      case TokenType::Times:
      {
        NodeTimes* node = new NodeTimes (nargs, tk.vfloc);
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
        NodeTimes* node = new NodeTimes (nargs, tk.vfloc);
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
        NodePower* node = new NodePower (tk.vfloc);
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
        NodeEquals* node = new NodeEquals (tk.vfloc);
        // put the argument from the stack
        node->addArgument (1, ts.back());
        ts.pop_back ();
        node->addArgument (0, ts.back());
        ts.pop_back ();
        ts.push_back (node);
      }
      break;
      case TokenType::Define:
      {
        NodeEquals* node = new NodeEquals (true, tk.vfloc);
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
        NodeSemicolon* node = new NodeSemicolon (tk.vfloc);
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
        PE_MESSAGE1(tk.vfloc, "Unexpected operator.");
        break;
    }
  } else PE_MESSAGE1(tk.vfloc, "Too short stack.");
}

Node* toPostfix (const std::vector<Token>& token_stream, const std::map<TokenType, TokenPrecedence>& token_table)
{
  // this is the current stack of operations
  // what can get onto this stack
  // - functions
  // - operators, until a higher presedence comes along
  // - opening parenthesis
  // everyting else goes straight onto the list for the tree
  std::list<Token> ops;
  // the size of this is the depth of function calls
  // the value at the given level is the number of function arguments
  std::list<size_t> funargs;
  // the tree consists sequence of [list of arguments followed by the {operator or function}]
  std::list<Node*> treestack;
  try
  {
    size_t k = 0;
    while (k < token_stream.size())
    {
      if (token_stream[k].type == TokenType::Symbol)
      {
        // if there is a funtion call and if there is no argument yet, 
        // then this will be the first argument
        if (!funargs.empty()) if (funargs.back() == 0) funargs.back() = 1;
        // bind items on the stack ...
        opToTree (treestack, token_stream[k], 0);
      }
      else if (token_stream[k].type == TokenType::Function)
      {
        // if there is a funtion call and if there is no argument yet, 
        // then this will be the first argument
        if (!funargs.empty()) if (funargs.back() == 0) funargs.back() = 1;
        ops.push_back (token_stream[k]);
      }
      else if (token_stream[k].type == TokenType::Comma)
      {
        if (!funargs.empty()) 
        {
          if (funargs.back() == 0)
          {
            // if there is a funtion call and if there is no argument yet, 
            // then it means we have an empty argument
            PE_MESSAGE1(token_stream[k].vfloc, "Inappropriate use of comma: first function argument is missing.\n");
          }
          else 
          {
            // we add an argument, if we already had one
            funargs.back() += 1;
          }
        } 
        else
        {
          // if this is not a function
          PE_MESSAGE1(token_stream[k].vfloc, "Inappropriate use of comma: it is not a function.\n");
        }
        // the comma seperates expression, so things can be moved into the tree
        // -> until there is a paranthesis (cannot be another comma, because we have emptied that before
        while ( (!ops.empty())&&(ops.back().type != TokenType::LeftParen) )
        {
          // there cannot be a function, but everyting else is allowed
          // if there is a function then we should have detected it???
          PE_ERROR_X1(ops.back().type != TokenType::Function, ops.back().vfloc, "Internal parser error: unexpected function.");
          opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
          ops.pop_back();
        }
        // it is wrong to have nothing left
        if (ops.empty()) PE_MESSAGE1(token_stream[k].vfloc, "Misplaced comma or missing \')\'\n");
      }
      // what happens with an operator
      else if (token_stream[k].isOp())
      {
        if (!funargs.empty()) if (funargs.back() == 0) funargs.back() = 1;
        while (!ops.empty())
        {
          // - it must be an operator
          // - if left associative, it must have lower or equal precedence then the current operator
          // - if not left associative it must have strictly lower then the current operator 
          if ( (ops.back().isOp())&&
              (((token_table.at(token_stream[k].type).leftAssoc)&&
                (token_table.at(token_stream[k].type).precedence >= token_table.at(ops.back().type).precedence))||
                (token_table.at(token_stream[k].type).precedence > token_table.at(ops.back().type).precedence)) )
          {
            // bind items on the stack ...
            PE_ERROR_X1(ops.back().type != TokenType::Function, ops.back().vfloc, "Internal parser error: unexpected function.");
            opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
            ops.pop_back();
          } else break;
        }
  //      std::cout <<"op(" << token_stream[k].name << ")=" << funargs.back () << " ";
        ops.push_back(token_stream[k]);
      }
      else if (token_stream[k].type == TokenType::LeftParen)
      {
  //      std::cout <<"Left : " << funargs.back () << " ";
        // any left paren starts to count the function arguments
        funargs.push_back(0);
        ops.push_back(token_stream[k]);
      }
      // new sub expression is opened
      else if (token_stream[k].type == TokenType::RightParen)
      {
        // this assumes that we are in a function context
        PE_ERROR_X1(!funargs.empty(), token_stream[k].vfloc, "Internal parser error: not in a function context");
        size_t currentArgs = funargs.back ();
  //      std::cout <<"Right : " << funargs.back () << " ";
        size_t it = 0;
        while (!ops.empty())
        {
          if (ops.back().type != TokenType::LeftParen)
          {
            // push all operations (ops) until a left parenthesis comes
            if (ops.back().type == TokenType::Function) opToTree (treestack, ops.back(), currentArgs);
            else opToTree (treestack, ops.back(), token_table.at(ops.back().type).nargs);
            ops.pop_back();
            it++;
          } else
          {
            // just discard the left parenthesis
            funargs.pop_back();
            ops.pop_back();
            break;
          }
        }
        // if nothing was put onto the tree stack and we are in function context then there is one argument
        if ((!funargs.empty()) && (it != 0)) if (funargs.back() == 0) funargs.back() = 1;
        if (ops.empty())
        {
          // there was noting on the stack
          std::ostringstream msg;
          msg << "Mismatched prenthesis (2)\nBEGIN@ ";
          for (size_t p=0; p<k+1; p++) msg << token_stream[p].name;
          msg << " @ ";
          for (size_t p=k+1; p<token_stream.size(); p++) msg << token_stream[p].name;
          msg << " @END\n ";
          PE_MESSAGE1(token_stream[k].vfloc, msg.str().c_str());
        }
        if (!ops.empty()) if (ops.back().type == TokenType::Function)
        {
          // bind items on the stack ...
          opToTree (treestack, ops.back(), currentArgs);
          ops.pop_back();
        }
      } else
      {
        PE_MESSAGE3(ops.back().vfloc, "Syntax error: unexpected token.", ops.back().name, "\n");
      }
      k++;
    }
    // no more tokens
    while (!ops.empty())
    {
      if (ops.back().type == TokenType::LeftParen) PE_MESSAGE1(ops.back().vfloc, "Mismatched prenthesis (3).");
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
      msg << "Syntax error: could not parse " << treestack.size() << "\n";
      for (auto k=treestack.begin(); k != treestack.end(); ++k)
      {
        msg << "@@OP =\n\t";
        if (*k) { (*k)->print (msg); msg << "\n"; }
      }
      for (auto k=token_stream.begin(); k != token_stream.end(); ++k)
      {
        msg << k->name << "{" << static_cast<int>(k->type) << "}";
      }
      P_MESSAGE1(msg.str().c_str());
    }
  }
  catch (KNException& ex)
  {
    while (!treestack.empty())
    {
      Node* it = treestack.back ();
      it->deleteTree();
      delete it;
      treestack.pop_back ();
    }
    throw ex;
  }
  if (treestack.empty()) return nullptr;
  else return treestack.front ();
}

// should use a map for par_name and par_value

// needs to handle period(), time(), PAR_PERIOD
void splitExpression (std::string& sysName,
  std::vector<NodeSymbol*>& var_name,
  std::vector<Node*>& var_dot,
  std::vector<Node*>& var_init,
  std::vector<Node*>& var_mass,
  std::vector<NodeSymbol*>& par_name,
  std::vector<Node*>& par_value,
  std::vector<NodeSymbol*>& time,
  std::vector<Node*>& macro,
  Node* root)
{
  // set this up for the period
  par_name.push_back(nullptr);
  par_value.push_back(nullptr);
  NodeSemicolon* parent = dynamic_cast<NodeSemicolon*>(root);
  if (parent != nullptr)
  {
    Node* cp = parent -> copy ();
    std::list<Node*> lst;
    static_cast<NodeSemicolon*>(cp)->toList(lst);
    try {
      for (auto it : lst)
      {
        NodeEquals* child = dynamic_cast<NodeEquals*>(it);
        if (child != nullptr)
        {
          const Node* left = child->getChild (0);
          const Node* right = child->getChild (1);

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
                  var_init.push_back (nullptr);
                  var_mass.push_back (nullptr);
                }
              } else PE_MESSAGE1(fun->vfloc, "Dot: wrong number of arguments\n");
            } else if (fun->getName() == "init")
            {
              // starting solution
              if (fun->children() == 1)
              {
                const Node* init = fun->getChild(0);
                if (init->type == TokenType::Symbol)
                {
                  const NodeSymbol* sym = static_cast<const NodeSymbol*>(init);
                  bool found = false;
                  // find the corresponding state variable
                  for (size_t k=0; k<var_name.size(); k++)
                  {
                    if (var_name[k]->getName() == sym->getName())
                    {
                      var_init[k] = right->copy();
                      found = true;
                      break;
                    }
                  }
                  PE_ERROR_X3(found == true, sym->vfloc, "State variable '", sym->getName(), "' is not defined.");
                } else PE_MESSAGE1(init->vfloc, "Init: state variable is not a symbol.");
              } else PE_MESSAGE1(fun->vfloc, "Init: wrong number of arguments\n");
            } else if (fun->getName() == "mass")
            {
              // starting solution
              if (fun->children() == 1)
              {
                const Node* init = fun->getChild(0);
                if (init->type == TokenType::Symbol)
                {
                  const NodeSymbol* sym = static_cast<const NodeSymbol*>(init);
                  bool found = false;
                  // find the corresponding state variable
                  for (size_t k=0; k<var_name.size(); k++)
                  {
                    if (var_name[k]->getName() == sym->getName())
                    {
                      var_mass[k] = right->copy();
                      found = true;
                      break;
                    }
                  }
                  PE_ERROR_X3(found == true, sym->vfloc, "State variable '", sym->getName(), "' is not defined");
                } else PE_MESSAGE1(init->vfloc, "Mass: state variable is not a symbol");
              } else PE_MESSAGE1(fun->vfloc, "Mass: wrong number of arguments");
            } else if (fun->getName() == "period")
            {
              if (par_name.size() > 0) PE_ERROR_X1(par_name[0] == nullptr, fun->vfloc, "period is already defined");
              PE_ERROR_X1(fun->children() == 1, fun->vfloc, "period: wrong number of arguments");
              PE_ERROR_X1(fun->getChild(0)->type == TokenType::Symbol, fun->getChild(0)->vfloc, "Not a symbol");
              par_name[0] = static_cast<NodeSymbol*>(fun->getChild(0)->copy());
              par_value[0] = right->copy();
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
              } else PE_MESSAGE1(fun->vfloc, "Par: wrong number of arguments\n");
            } else if (fun->getName() == "time")
            {
              PE_ERROR_X1(fun->children() == 0, fun->vfloc, "Too many arguments");
              PE_ERROR_X1(right->type == TokenType::Symbol, fun->vfloc, "Not a symbol");
              time.push_back (static_cast<NodeSymbol*>(right->copy()));
            } else if (fun->getName() == "vfname")
            {
              PE_ERROR_X1(fun->children() == 0, fun->vfloc, "Too many arguments");
              PE_ERROR_X1(right->type == TokenType::Symbol, fun->vfloc, "Not a symbol");
              sysName = static_cast<const NodeSymbol*>(right)->getName ();
            } else
            {
              // a simple macro
              macro.push_back(child->copy());
            }
          } else if (left->type == TokenType::Symbol)
          {
            macro.push_back(child->copy());
          } else
          {
            left->print (std::cout); std::cout << " == "; right->print (std::cout); std::cout << "\n";
            std::cout.flush();
            PE_MESSAGE1(child->vfloc, "Not a definition\n");
          }
        } else
        {
          PE_MESSAGE1(it->vfloc, "Not a definition\n");
        }
      }
    }
    catch (KNException& ex)
    {
      for (auto it : lst) { it->deleteTree (); delete it; }
      cp -> deleteTree ();
      delete cp;
      throw ex;
    }
    for (auto it : lst) { it->deleteTree (); delete it; }
    cp -> deleteTree ();
    delete cp;
  } else
  {
    PE_MESSAGE1(root->vfloc, "Not a semicolon, cannot split expression");
  }
  if (par_name[0] == nullptr)
  {
    par_name[0] = new NodeSymbol ("T", 0);
    par_value[0] = new NodeNumber (1.0, 0);
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
    token_table[TokenType::Define] = { 10, false, 2 };
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
  std::vector<std::string>& exprName,
  std::vector<Expression>& exprFormula,
  std::vector<Expression>& varInit,
  std::vector<double>& varMass,
  std::vector<Expression>& delayExpr,
  std::vector<std::string>& parName,
  std::vector<double>& parInit)
{
  std::vector<NodeSymbol*> var_name;
  std::vector<Node*> var_dot;
  std::vector<Node*> var_init;
  std::vector<Node*> var_mass;
  std::vector<NodeSymbol*> par_name;
  std::vector<Node*> par_value;
  std::vector<NodeSymbol*> time;
  std::vector<Node*> all_def;
  std::vector<Node*> macro;

  // Adding Pi
  NodeEquals* pi_eq = new NodeEquals (0);
  pi_eq -> addArgument (0, new NodeSymbol("Pi", 0));
  pi_eq -> addArgument (1, new NodeNumber(M_PI, 0));
  macro.push_back (pi_eq);

  try {
    splitExpression (sysName, var_name, var_dot, var_init, var_mass, par_name, par_value, time, all_def, root);
  }
  catch (KNException& ex)
  {
    // cleaning up the variables
    for (auto it : var_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_dot) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_init) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_mass) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : par_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : par_value) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : time) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : macro) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : all_def) { if (it != nullptr) { it->deleteTree(); delete it; } }
//     std::cout << "From knutSplit - ";
    throw ex;
  }
  // need to further split macros into replacements and definitions
  // transfers ownership from all_ref
  std::vector<NodeSymbol*> expr_name;
  std::vector<Node*> expr_formula;
  for (auto it : all_def)
  {
    if (it->type == TokenType::Define)
    {
      NodeSymbol* sym = dynamic_cast<NodeSymbol*>(it->getChild(0));
      if (sym != nullptr)
      {
        expr_name.push_back (sym);
        expr_formula.push_back (it->getChild(1));
        delete it;
      }
      else
      {
        std::cout << "Invalid expression: ";
        it->print (std::cout);
        std::cout << std::endl;
        it->deleteTree ();
        delete it;
      }
    } else
    {
      macro.push_back (it);
    }
  }
  all_def.resize (0);

  // if no independent variable or period is set we set the default
  if (time.empty()) time.push_back (new NodeSymbol ("t", 0));
  NodeVar* time_var = new NodeVar (0, 0);

  // adding internal parameters
  par_name.push_back (new NodeSymbol ("PAR_ANGLE", 0));
  par_value.push_back (new NodeNumber (0.0, 0));
  par_name.push_back (new NodeSymbol ("PAR_PERIOD", 0));
  par_value.push_back (new NodeNumber (1.0, 0));
  par_name.push_back (new NodeSymbol ("PAR_ROT_NUM", 0));
  par_value.push_back (new NodeNumber (0.0, 0));

  std::vector<NodeVar*> var_idx;
  std::vector<NodePar*> par_idx;
  std::vector<NodeExpr*> expr_idx;

  // dummy variables and parameters
  std::vector<double> var_numval;

  varName.resize (var_name.size());
  for (size_t k = 0; k < var_name.size(); k++)
  {
    var_idx.push_back (new NodeVar (1 + k, 0)); // idx == 0 is the time
    varName[k] = var_name[k]->getName ();
  }

  parName.resize (par_name.size());
  for (size_t k = 0; k < par_name.size(); k++)
  {
    par_idx.push_back (new NodePar (k, 0)); // idx == 0 is the period
    parName[k] = par_name[k]->getName ();
  }

  exprName.resize (expr_name.size());
  for (size_t k = 0; k < expr_name.size(); k++)
  {
    expr_idx.push_back (new NodeExpr (k, 0));
    exprName[k] = expr_name[k]->getName ();
  }

  // Applying the macros to both parameters and variables
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
    for (size_t q = 0; q < macro.size (); q++)
    {
      par_value[p] -> replaceSymbol (*macro[q], &par_value[p]);
    }
  }

  for (size_t p = 0; p < macro.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      if (p != q) macro[p] -> replaceSymbol (*macro[q], &macro[p]);
    }
  }

  // var
  for (size_t p = 0; p < var_dot.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      var_dot[p] -> replaceSymbol (*macro[q], &var_dot[p]);
    }
  }

  // expr
  for (size_t p = 0; p < expr_formula.size(); p++)
  {
    for (size_t q = 0; q < macro.size(); q++)
    {
      expr_formula[p] -> replaceSymbol (*macro[q], &expr_formula[p]);
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

  for (size_t p = 0; p < var_mass.size(); p++)
  {
    if (var_mass[p] != nullptr)
    {
      for (size_t q = 0; q < macro.size(); q++)
      {
        var_mass[p] -> replaceSymbol (*macro[q], &var_mass[p]);
      }
      for (size_t q = 0; q < par_name.size(); q++)
      {
        var_mass[p] -> replaceSymbol (*par_name[q], *par_value[q], &var_mass[p]);
      }
      var_mass[p] -> replaceSymbol (*time[0], *time_var, &var_mass[p]);
    }
  }

  // fill up the parameter values
  parInit.resize (par_name.size());
  for (size_t p = 0; p < par_name.size(); p++)
  {
    par_value[p] -> evaluateWithPar (&parInit[p], parInit);
  }

  // fill up the starting point
  varInit.resize (var_init.size());
  for (size_t p = 0; p < var_init.size(); p++)
  {
    varInit[p].fromNode (var_init[p]);
    var_init[p] = nullptr;
  }
  var_init.resize (0);

  // fill up the mass
  varMass.resize (var_mass.size ());
  for (size_t p = 0; p < var_mass.size (); p++)
  {
    if (var_mass[p] != nullptr)
    {
      var_mass[p] -> evaluateWithPar (&varMass[p], parInit);
    } else varMass[p] = 1.0;
  }

  // START :: HANDLING THE DELAYS
  // do not delete any of these since they are owned by the tree!
  std::vector<const NodeFunction*> delays;
  for (size_t p = 0; p < var_dot.size(); p++)
  {
    var_dot[p] -> findFunction (delays, "delay");
  }
  for (size_t p = 0; p < expr_formula.size(); p++)
  {
    expr_formula[p] -> findFunction (delays, "delay");
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
      NodeVar var_del(1 + q + var_name.size()*(delay_id[p]+1), 0);
      fun->replaceSymbol (*var_name[q], var_del, nullptr);
    }
  }

  var_numval.resize (var_name.size()*(1+unique_delays.size()));
  // fill up the delays
  delayExpr.resize (unique_delays.size());
  Node* delay = nullptr;
  try
  {
    for (size_t r = 0; r < unique_delays.size(); r++)
    {
      delay = unique_delays[r]->getChild (1)->copy();
      // replace parameters
      for (size_t q = 0; q < par_name.size(); q++)
      {
        delay -> replaceSymbol (*par_name[q], *par_idx[q], &delay);
      }
      // replace time
      delay -> replaceSymbol (*time[0], *time_var, &delay);
      // take over ownership
      delayExpr[r].fromNode (delay);
      delay = nullptr;
//      // evaluate to find undefined symbols
//      delayExpr[r].test ();
    }
  }
  catch (KNException& ex)
  {
    if (delay != nullptr) { delay->deleteTree(); delete delay; }
    delayExpr.resize (0);
//     std::cout << "From knutSplit(delays) - ";
    throw ex;
  }

  // render delay into identity
  NodeEquals* del_rem = new NodeEquals (0);
  NodeFunction* del_id = new NodeFunction ("delay", 2, 0);
  del_id->addArgument (0, new NodeSymbol ("a", 0));
  del_id->addArgument (1, new NodeSymbol ("b", 0));
  del_rem->addArgument (0, del_id);
  del_rem->addArgument (1, new NodeSymbol ("a", 0));

  // replacing the variables
  // in dot
  for (size_t p = 0; p < var_dot.size(); p++)
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
    for (size_t q = 0; q < expr_name.size(); q++)
    {
      var_dot[p] -> replaceSymbol (*expr_name[q], *expr_idx[q], &var_dot[p]);
    }
  }
  // in expr
  for (size_t p = 0; p < expr_formula.size(); p++)
  {
    expr_formula[p] -> replaceSymbol (*del_rem, &expr_formula[p]);
    expr_formula[p] -> replaceSymbol (*time[0], *time_var, &expr_formula[p]);
    for (size_t q = 0; q < var_name.size(); q++)
    {
      expr_formula[p] -> replaceSymbol (*var_name[q], *var_idx[q], &expr_formula[p]);
    }
    for (size_t q = 0; q < par_name.size(); q++)
    {
      expr_formula[p] -> replaceSymbol (*par_name[q], *par_idx[q], &expr_formula[p]);
    }
    // only until the current definition so that no recursion is possible
    for (size_t q = 0; q < p; q++)
    {
      expr_formula[p] -> replaceSymbol (*expr_name[q], *expr_idx[q], &expr_formula[p]);
    }
  }
  del_rem -> deleteTree();
  delete del_rem;
  delete time_var;
  // END :: HANDLING THE DELAYS

  // fill up the RHS
  try
  {
    varDotExpr.resize (var_dot.size());
    for (size_t q = 0; q < var_dot.size(); q++)
    {
      // transfer ownership
      varDotExpr[q].fromNode (var_dot[q]);
      var_dot[q] = nullptr;
//      // evaluate to find undefined symbols
//      varDotExpr[q].test ();
    }
    var_dot.resize (0);
    // fill up the expressions
    exprFormula.resize (expr_formula.size());
    for (size_t q = 0; q < expr_formula.size(); q++)
    {
      // transfer ownership
      exprFormula[q].fromNode (expr_formula[q]);
      expr_formula[q] = nullptr;
//      // evaluate to find undefined symbols
//      exprFormula[q].test ();
    }
    expr_formula.resize (0);
  }
  catch (KNException& ex)
  {
    varDotExpr.resize (0);
    exprFormula.resize (0);
//     std::cout << "From knutSplit(varDotExpr, exprFormula) - ";
    for (auto it : var_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_dot) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_init) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_mass) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : par_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : par_value) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : time) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : macro) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : expr_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : expr_formula) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : var_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : par_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
    for (auto it : expr_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
    throw ex;
  }
  // cleaning up the variables
  for (auto it : var_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : var_dot) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : var_init) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : var_mass) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : par_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : par_value) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : time) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : macro) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : expr_name) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : expr_formula) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : var_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : par_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
  for (auto it : expr_idx) { if (it != nullptr) { it->deleteTree(); delete it; } }
}

#include <mxml.h>

bool isValidName(const std::string& s)
{
  const std::string lower("abcdefghijklmnopqrstuvwxyz");
  const std::string upper("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  const std::string letters = lower + upper + "_(),";
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

// This is a compatibility function, it is depreceated
// converts xmlstring into an expression to be parsed.
// returns true if it was a real XML file, otherwise, false
bool Expression::fromXML (std::string& oexpr, const std::string& xmlstring)
{
  std::ostringstream out;
  mxml_node_t *tree;
  mxml_node_t *node;
//   bool bad_attr;

  tree = mxmlLoadString (nullptr, xmlstring.c_str (), MXML_NO_CALLBACK);
  if (tree == nullptr)
  {
    // Not an XML string
    return false;
  }

  node = mxmlFindElement(tree, tree, "VectorField", nullptr, nullptr, MXML_DESCEND);
  if (node == nullptr)
  {
    mxmlDelete (tree);
//     P_MESSAGE1("Error: No VectorField element found in XML defintion.\n");
    return false;
  }
  else
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "IndependentVariable" && attr != "Description")
      {
        mxmlDelete (tree);
        P_MESSAGE2("Error: The VectorField element has an unknown attribute: ", attr);
      }
    }
    const char *attr = mxmlElementGetAttr(node, "Name");
    if (attr == nullptr)
    {
      mxmlDelete (tree);
      P_MESSAGE1("Error: The VectorField element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The VectorField Name \"", attr, "\" is not valid.");
      }
      out << "vfname()=" << attr << ";\n";
    }
    attr = mxmlElementGetAttr(node, "Description");
    if (attr != nullptr)
    {
      out << "vfdescription()=" << attr << ";\n";
    }
    attr = mxmlElementGetAttr(node, "IndependentVariable");
    if (attr == nullptr)
    {
      out << "time()=t;\n";
    } else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The VectorField IndependentVariable \"", attr, "\" is not valid.");
      }
      out << "time()=" << attr << ";\n";
    }
  }

  //
  // Get the constants
  //
  for (node = mxmlFindElement(tree, tree, "Constant", nullptr, nullptr, MXML_DESCEND);
       node != nullptr;
       node = mxmlFindElement(node, tree, "Constant", nullptr, nullptr, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Value" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: A Constant element has an unknown attribute: ", attr, "Valid Constant attributes are: Name, Value, Description, Latex.");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == nullptr)
    {
      mxmlDelete (tree);
      P_MESSAGE1("Error: A Constant element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The Constant Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << name << "=";
      attr = mxmlElementGetAttr(node, "Value");
      if (attr == nullptr)
      {
        mxmlDelete (tree);
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
  for (node = mxmlFindElement(tree, tree, "Parameter", nullptr, nullptr, MXML_DESCEND);
       node != nullptr;
       node = mxmlFindElement(node, tree, "Parameter", nullptr, nullptr, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "DefaultValue" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: A Parameter element has an unknown attribute: ", attr, "Valid Parameter attributes are: Name, DefaultValue, Description, Latex.");
//         bad_attr = true;
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == nullptr)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: A Parameter element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The Parameter Name \"", attr, "\" is not valid.\n");
      }
      std::string name(attr);
      std::string descr;
//       Parameter *p = new Parameter(name);
      attr = mxmlElementGetAttr(node, "Description");
      if (attr != nullptr)
      {
        descr = attr;
        std::transform(descr.begin(), descr.end(), descr.begin(), ::tolower);
      }
      if (descr == "period")
        out << "period(" << name << ")=";
      else
        out << "par(" << name << ")=";

      attr = mxmlElementGetAttr(node, "DefaultValue");
      if (attr != nullptr)
      {
        out << attr << ";\n";
      }
    }
  }

  //
  // Get the auxiliary expressions
  //
  for (node = mxmlFindElement(tree, tree, "Expression", nullptr, nullptr, MXML_DESCEND);
       node != nullptr;
       node = mxmlFindElement(node, tree, "Expression", nullptr, nullptr, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Formula" && attr != "Description" && attr != "Latex")
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: An Expression element has an unknown attribute: ", attr, "Valid Expression attributes are: Name, Formula, Description, Latex");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == nullptr)
    {
      mxmlDelete(tree);
      P_MESSAGE1("Error: An Expression element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The Expression Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << name << "=";
      attr = mxmlElementGetAttr(node, "Formula");
      if (attr == nullptr)
      {
        mxmlDelete (tree);
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
  for (node = mxmlFindElement(tree, tree, "StateVariable", nullptr, nullptr, MXML_DESCEND);
       node != nullptr;
       node = mxmlFindElement(node, tree, "StateVariable", nullptr, nullptr, MXML_DESCEND))
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
    if (attr == nullptr)
    {
      mxmlDelete (tree);
      P_MESSAGE1("Error: A StateVariable element has no Name attribute.");
    }
    else
    {
      if (!isValidName(attr))
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The StateVariable Name \"", attr, "\" is not valid.");
      }
      std::string name(attr);
      out << "dot(" << name << ")=";
      attr = mxmlElementGetAttr(node, "Formula");
      if (attr == nullptr)
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: The StateVariable with Name=\"", name, "\" has no Formula attribute.");
      }
      else
      {
        out << attr << ";\n";
      }
      attr = mxmlElementGetAttr(node, "DefaultInitialCondition");
      if (attr != nullptr)
      {
        out << "init(" << name << ")=" << attr << ";\n";
      }
      attr = mxmlElementGetAttr(node, "Mass");
      if (attr != nullptr)
      {
        out << "mass(" << name << ")=" << attr << ";\n";
      }
    }
  }

  //
  // Get the functions
  //
  for (node = mxmlFindElement(tree, tree, "Function", nullptr, nullptr, MXML_DESCEND);
       node != nullptr;
       node = mxmlFindElement(node, tree, "Function", nullptr, nullptr, MXML_DESCEND))
  {
    for (int i = 0; i < node->value.element.num_attrs; ++i)
    {
      std::string attr = node->value.element.attrs[i].name;
      if (attr != "Name" && attr != "Formula" && attr != "Description")
      {
        mxmlDelete (tree);
        P_MESSAGE3("Error: A Function element has an unknown attribute: ", attr, "Valid Function attributes are: Name, Formula, Description.");
      }
    }
    const char *attr;
    attr = mxmlElementGetAttr(node, "Name");
    if (attr == nullptr)
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
      if (attr == nullptr)
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

  mxmlDelete (tree);
  oexpr = out.str ();
  return true;
}

} // namespace
