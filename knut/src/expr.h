// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 - 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef EXPR_H
#define EXPR_H

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <list>
#include <functional>
#include <ostream>
#include <iostream>
#include "knerror.h"

#define KNUT_STACK_INCREMENT 5

namespace ExpTree
{

enum class TokenType : int
{
  Invalid=0,
  Function,
  UnaryPlus,
  UnaryMinus,
  Plus,
  Minus,
  Times,
  Divide,
  Power,
  LeftParen,
  RightParen,
  Comma,
  Equals,
  Define,
  Semicolon,
  Symbol,
  Number,
  Var,
  Par,
  Expr
};

// the pointer is owned somewhere else
class Value
{
public:
  Value() : data(nullptr), skip(1) { }
  ~Value() { }
  Value (const Value& value) { data = value.data; skip = value.skip; }
  double* data;
  size_t skip;
};

class ConstValue
{
public:
  ConstValue() : data(nullptr), skip(1) { }
  ~ConstValue() { }
  ConstValue (const ConstValue& value) { data = value.data; skip = value.skip; }
  ConstValue (const Value& value) { data = value.data; skip = value.skip; }
  const double* data;
  size_t skip;
};

class ValueStack
{
public:
  ValueStack ();
  ValueStack (double *res, size_t w) : width(w), stack(1)
  {
    stack[0].data = res;
    stack[0].skip = 1;
  }
  ~ValueStack ()
  {
    for (size_t k = 1; k < stack.size(); k++) delete[] stack [k].data;
  }
  void setResult (double *res, size_t skip = 1)
  {
    stack[0].data = res;
    stack[0].skip = skip;
  }
  Value& operator[](size_t p) { return stack[p]; }
  void resizeDepth (size_t depth)
  {
    const size_t from = stack.size();
    if (depth > from) stack.resize(depth);
    for (size_t k = from; k < depth; k++)
    {
      stack [k].data = new double[width];
      stack [k].skip = 1;
    }
  }
  void resizeWidth (size_t w)
  {
    if (w > width)
    {
      for (size_t k = 1; k < stack.size(); k++)
      {
        delete[] stack [k].data;
        stack [k].data = new double[w];
        stack [k].skip = 1;
      }
      width = w;
    }
  }
  void resize (size_t depth, size_t w)
  {
    const size_t from = stack.size();
    const size_t nw = (w > width) ? w : width;
    if (depth > from) stack.resize(depth);
    for (size_t k = from; k < depth; k++)
    {
      stack [k].data = new double[nw];
      stack [k].skip = 1;
    }
    if (w > width)
    {
      for (size_t k = 1; k < from; k++)
      {
        delete[] stack [k].data;
        stack [k].data = new double[w];
        stack [k].skip = 1;
      }
      width = w;
    }
  }
  size_t size() const { return stack.size(); }
private:
  size_t width;
  std::vector<Value> stack;
};

// forward declaration
class NodeSymbol;
class NodeFunction;

class Node
{
public:
  Node (TokenType tp, const size_t loc) : type(tp), vfloc(loc)
  {
#ifdef DEBUG
    instances.push_front(this);
    ist = instances.begin();
#endif
  }
  virtual ~Node()
  {
#ifdef DEBUG
    instances.erase (ist);
#endif
  }
  virtual void deleteTree () = 0;
  virtual Node* copy () const = 0;
  virtual size_t children () const = 0;
  virtual const Node* getChild (size_t n) const = 0;
  virtual Node* getChild (size_t n) = 0;
  virtual size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const = 0;
//  virtual void addArgument (size_t pos, Node* arg) = 0;
  virtual void print (std::ostream& out) const = 0;
  virtual bool operator != (const Node& node) const = 0;
  bool operator == (const Node& node) const { return !((*this) != node); }
  virtual bool operator < (const Node& node) const = 0;
  virtual void replaceSymbol(const Node& sym, const Node& node, Node** parent) = 0;
  virtual void replaceSymbol(const Node& node, Node** parent) = 0;
  virtual void findFunction(std::vector<const NodeFunction*>& lst, const std::string& name) const = 0;
  virtual Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const = 0;
  // optimize the node
  // return: this (pointer to the node itself : nothing to do
  // return: anything else: delete the node and replace with the new one
  virtual Node* optimize (const std::function<bool(const Node*)>& zero) = 0;
  virtual void find (TokenType tp, const std::function<void(const Node*)>& found) const;
  void evaluateWithPar (double *res, const std::vector<double>& parInit);
  const TokenType type;
  const size_t vfloc; // location in the .vf file
  static const std::function<bool(const Node*)> alwaysFalse;
#ifdef DEBUG
  static std::list<const Node *> instances;
  std::list<const Node *>::iterator ist;
#endif
  static inline Node* node_optimize(Node* nd, const std::function<bool(const Node*)>& zero)
  {
    Node* ret = nd -> optimize (zero);
    if (ret != nd)
    {
      delete nd;
      return ret;
    }
    return nd;
  }
};

static inline Node* node_optimize(Node* nd, const std::function<bool(const Node*)>& zero)
{
  Node* ret = nd -> optimize (zero);
  if (ret != nd)
  {
    delete nd;
    return ret;
  }
  return nd;
}

class NodeNumber : public Node
{
public:
  NodeNumber (double val, const size_t loc) : Node(TokenType::Number, loc), value(val) {}
  void deleteTree () override {}
  Node* copy () const override { return new NodeNumber (value, vfloc); }
  size_t children () const override { return 0; }
  const Node* getChild (size_t n) const override { return nullptr; }
  Node* getChild (size_t n) override { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override { if (value > 0) out << value; else out << '(' << value << ')'; }
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override {}
  void replaceSymbol(const Node& node, Node** parent) override {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const override {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override { return this; }
  double getValue() const { return value; }
private:
  double value;
};

class NodeSymbol : public Node
{
public:
  NodeSymbol (const std::string& nm, const size_t loc) : Node(TokenType::Symbol, loc), name(nm) {}
  void deleteTree () override {}
  Node* copy () const override { return new NodeSymbol (name, vfloc); }
  size_t children () const override { return 0; }
  const Node* getChild (size_t n) const override { return nullptr; }
  Node* getChild (size_t n) override { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override { out << name; }
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol (const Node& node, Node** parent) override;
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const override {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override { return this; }
  const std::string& getName () const { return name; }
private:
  std::string name;
};

class NodeVar : public Node
{
public:
  NodeVar (size_t id, const size_t loc) : Node(TokenType::Var, loc), idx(id) {}
  ~NodeVar() override {}
  void deleteTree () override {}
  Node* copy () const override { return new NodeVar (idx, vfloc); }
  size_t children () const override {return 0; }
  const Node* getChild (size_t n) const override { return nullptr; }
  Node* getChild (size_t n) override { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override { out << "V["<<idx<<"]"; }
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol (const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol (const Node& node, Node** parent) override {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const override {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override { return this; }
  size_t getIdx () const { return idx; }
private:
  size_t idx;
};

class NodePar : public Node
{
public:
  NodePar (size_t id, const size_t loc) : Node(TokenType::Par, loc), idx(id) {}
  ~NodePar() override {}
  void deleteTree () override {}
  Node* copy () const override { return new NodePar (idx, vfloc); }
  size_t children () const override { return 0; }
  const Node* getChild (size_t n) const override { return nullptr; }
  Node* getChild (size_t n) override { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override { out << "P["<<idx<<"]"; }
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const override {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override { return this; }
  size_t getIdx () const { return idx; }
private:
  size_t idx;
};

class NodeExpr : public Node
{
public:
  NodeExpr (size_t id, const size_t loc) : Node(TokenType::Expr, loc), idx(id), np(0), nv(0) { }
  NodeExpr (size_t idx_, size_t np_, size_t nv_, size_t dp_, size_t dv0_, size_t dv1_, const size_t loc)
    : Node(TokenType::Expr, loc), idx(idx_), np(np_), nv(nv_), dp(dp_) { dv[0] = dv0_; dv[1] = dv1_; }

  ~NodeExpr() override { }
  void deleteTree () override {}
  NodeExpr* copy () const override
  {
    NodeExpr *ne = new NodeExpr (idx, vfloc);
    ne->np = np;
    ne->nv = nv;
    ne->dp = dp;
    ne->dv[0] = dv[0];
    ne->dv[1] = dv[1];
    return ne;
  }
  size_t children () const override { return 0; }
  const Node* getChild (size_t n) const override { return nullptr; }
  Node* getChild (size_t n) override { return nullptr; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override
  {
    out << "ex" << idx;
    if (np > 0) out << "_p" << dp;
    if (nv == 1) out << "_v" << dv[0];
    if (nv == 2)
    {
      if (dv[0] < dv[1]) out << "_v" << dv[0] << "_v" << dv[1];
      else out << "_v" << dv[1] << "_v" << dv[0];
    }
  }
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override {}
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& name) const override {}
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override
  {
    if (zero (this))
    {
      Node* tmp = new NodeNumber (0.0, vfloc);
      deleteTree ();
      return tmp;
    }
    return this;
  }
  size_t getIdx () const { return idx; }
  size_t getNp () const { return np; }
  size_t getNv () const { return nv; }
  size_t getDp () const { return dp; }
  size_t getDv0 () const { return dv[0]; }
  size_t getDv1 () const { return dv[1]; }
private:
  size_t idx;
  size_t np;
  size_t nv;
  size_t dp;
  size_t dv[2];
};

class NodeFunction : public Node
{
public:
  NodeFunction (const std::string& nm, size_t nargs, const size_t loc) : Node(TokenType::Function, loc), name(nm), args(nargs)
  {
    for (size_t k=1; k<nargs; k++) args[k] = nullptr;
  }
  void deleteTree() override
  {
    for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree(); delete args[k]; }
  }
  Node* copy () const override;
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override { return args[n]; }
  Node* getChild (size_t n) override { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override;
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  // replaces the symbol with the node
  // parent points to the pointer of the current class
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override;
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override;
  // specific
  NodeFunction* specialize ();
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
  const std::string& getName() const { return name; }
protected:
  std::string name;
  std::vector<Node*> args;
};

class NodeAdd : public Node
{
public:
  NodeAdd (size_t nargs, const size_t loc) : Node(TokenType::Plus, loc), args(nargs), mul(nargs)
  {
    for (size_t k=1; k<nargs; k++) { args[k] = nullptr; mul[k] = 1.0; }
  }
  void deleteTree () override
  {
    for(size_t k=0; k<args.size(); k++)
      if (args[k] != nullptr) { args[k]->deleteTree (); delete args[k]; }
      else P_MESSAGE1("NULLPTR\n");
  }
  Node* copy () const override;
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override { return args[n]; }
  Node* getChild (size_t n) override { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override;
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override;
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  // specific
  void addArgument (size_t pos, Node* arg, double m) { args[pos] = arg; mul[pos] = m; }
  Node* optimize (const std::function<bool(const Node*)>& zero) override;
  void sort ();
private:
  std::vector<Node*> args;
  std::vector<double> mul;
};

class NodeTimes : public Node
{
public:
  NodeTimes (size_t nargs, const size_t loc) : Node(TokenType::Times, loc), args(nargs), smul(1.0), divide(nargs)
  {
    for (size_t k=1; k<nargs; k++) { args[k] = nullptr; divide[k] = false; }
  }
  void deleteTree () override
  {
    for(size_t k=0; k<args.size(); k++)
      if (args[k] != nullptr) { args[k]->deleteTree (); delete args[k]; }
      else P_MESSAGE1("NULLPTR\n");
  }
  Node* copy () const override;
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override { return args[n]; }
  Node* getChild (size_t n) override { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override;
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override;
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override;
  // needs to handle the divisions
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  // specific
  Node* optimize (const std::function<bool(const Node*)>& zero) override;
  void sort ();
  void setMul (double mul) { smul = mul; }
  double getMul () { return smul; }
  void addArgument (size_t pos, Node* arg, double m, bool div);
private:
  std::vector<Node*> args;
  double smul;
  std::vector<bool> divide;
};

class NodePower : public Node
{
public:
  NodePower (const size_t loc) : Node(TokenType::Power, loc), args(2) {}
  void deleteTree () override
  {
    for(size_t k=0; k<args.size(); k++) { args[k]->deleteTree (); delete args[k]; }
  }
  Node* copy () const override;
  size_t children () const override { return args.size(); }
  const Node* getChild (size_t n) const override { return args[n]; }
  Node* getChild (size_t n) override { return args[n]; }
  size_t evaluate (ValueStack& stack, size_t sp, const std::function<void(size_t, const Node*)>& fun, size_t len) const override;
  void print (std::ostream& out) const override;
  bool operator != (const Node& node) const override;
  bool operator < (const Node& node) const override;
  void replaceSymbol(const Node& sym, const Node& node, Node** parent) override;
  void replaceSymbol(const Node& node, Node** parent) override;
  void findFunction (std::vector<const NodeFunction*>& lst, const std::string& nm) const override;
  Node* derivative (const Node* var, const std::function<bool(const Node*)>& zero) const override;
  Node* optimize (const std::function<bool(const Node*)>& zero) override;
  // specific
  void addArgument (size_t pos, Node* arg) { args[pos] = arg; }
private:
  std::vector<Node*> args;
};

class Expression
{
public:
  Expression () : root (nullptr) { }
  Expression (const Expression& ex) // copy constructor necessary for containers
  {
    root = ex.root -> copy ();
  }
  Expression (const std::string& str);
  ~Expression ()
  {
    if (root != nullptr)
    {
      root -> deleteTree ();
      delete root;
    }
  }
  // copy operator
  const Expression& operator = (const Expression& exp); // disable copy operator
  // takes over the management of the Node rt
  void fromString (const std::string& str);
  void fromNode (Node* rt)
  {
    if (root != nullptr)
    {
      root -> deleteTree ();
      delete root;
    }
    root = rt;
  }
  const Node* node () const { return root; }
  void derivative (Expression& ex, const Node* nd,
                   const std::function<bool(const Node*)>& zero = Node::alwaysFalse) const
  {
    if(root != nullptr)
    {
      Node* deri = root->derivative (nd, zero);
      deri = node_optimize (deri, zero);
      ex.fromNode (deri);
    }
  }
  void optimize (const std::function<bool(const Node*)>& zero = Node::alwaysFalse)
  {
    if (root != nullptr) root = node_optimize (root, zero);
  }
  void print (std::ostream& out) const { if (root != nullptr) root -> print (out); }
  Node* copy () const { if(root != nullptr) return root -> copy (); else return nullptr; }
  void evaluate (ValueStack& stack, const std::function<void(size_t, const Node*)>& fun, size_t sp, size_t len) const
  {
    if (root != nullptr) root -> evaluate (stack, sp, fun, len);
  }
  bool isZero() const
  {
    if (root != nullptr)
    {
      const NodeNumber* num = dynamic_cast<NodeNumber*>(root);
      if (num != nullptr)
      {
        if (std::fabs(num -> getValue()) < 2*std::numeric_limits<double>::epsilon())
        {
          return true;
        } else return false;
      } else return false;
    } else return true;
  }
  size_t test ()
  {
    if (root != nullptr)
    {
//       double zero = 0.0;
      double res;
      ValueStack stk (&res, 1);
      auto fun = [&stk] (size_t sp, const Node*) { stk[sp].data[0] = 0.0; stk[sp].skip = 1; };
      root -> evaluate (stk, 0, fun, 1);
      return stk.size();
    }
    return 0;
  }
  void knutSplit (
    std::string& sysName,
    std::vector<std::string>& varName,
    std::vector<Expression>& varDotExpr,
    std::vector<std::string>& exprName,
    std::vector<Expression>& exprFormula,
    std::vector<Expression>& varInit,
    std::vector<double>& mass,
    std::vector<Expression>& delayExpr,
    std::vector<std::string>& parName,
    std::vector<double>& parInit);
  static bool fromXML (std::string& oexpr, const std::string& xmlstring);
private:
  // takes over ownership
//   Expression (Node* rt) : root(rt) { }
  Node* root;
};

} // namespace

#endif
