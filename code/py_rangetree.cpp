/*
 *
 * Program       : SRMCollider
 * Author        : Hannes Roest <roest@imsb.biol.ethz.ch>
 * Date          : 05.02.2011 
 *
 *
 * Copyright (C) 2011 - 2012 Hannes Roest
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
 *
 */

#include "rangetree.cpp"
#include "srmcolliderLib.cpp"
#include "py_srmcolliderLib.h"

// Expose to Python
using namespace python;
BOOST_PYTHON_MODULE(c_rangetree)
{

    def("create_tree", SRMCollider::SimpleRangetree::create_tree, 
 "Query the rangetree. Format is (x1,y1,x2,y2), returns all entries that are \n"
 "in the defined square defined by these four numbers. \n"
 "Returns a list with keys that were stored in the tree. \n"
 "\n"
 "\n"
 " Signature\n"
 "list query_tree(double a, double b, double c, double d)\n"
            "");
    def("query_tree", SRMCollider::SimpleRangetree::query_tree,
            
 "Create the rangetree that will be used throughout. This is essential. The \n"
 "rangetree will stay in place while this module is loaded. \n"
 "\n"
 "\n"
 " Signature\n"
 "void create_tree(tuple pepids) \n"
            "");

class_<SRMCollider::SimpleRangetree::Rangetree_Q1_RT,
  boost::shared_ptr<SRMCollider::SimpleRangetree::Rangetree_Q1_RT> >("Rangetree_Q1_RT",init<>())
        .def("create",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::create )
        .staticmethod("create")
        .def("new_rangetree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::new_rangetree)
        .def("create_tree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::create_tree)
        .def("query_tree",&SRMCollider::SimpleRangetree::Rangetree_Q1_RT::query_tree)
    ;

class_<SRMCollider::ExtendedRangetree::Rangetree_Q1_RT,
  boost::shared_ptr<SRMCollider::ExtendedRangetree::Rangetree_Q1_RT> >("ExtendedRangetree_Q1_RT",init<>())
        .def("create",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::create )
        .staticmethod("create")
        .def("new_rangetree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::new_rangetree)
        //.def("create_tree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::create_tree) // doesnt work any more since we overloaded create_tree
        .def("create_tree", static_cast<void (SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::*)(python::tuple)> (&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::create_tree) )
        .def("query_tree",&SRMCollider::ExtendedRangetree::Rangetree_Q1_RT::query_tree)
    ;
}

