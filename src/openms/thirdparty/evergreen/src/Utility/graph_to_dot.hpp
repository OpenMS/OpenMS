#ifndef _GRAPH_TO_DOT_HPP
#define _GRAPH_TO_DOT_HPP

#include <fstream>

template <typename VARIABLE_KEY>
void graph_to_dot(const InferenceGraph<VARIABLE_KEY> & ig, std::ostream & os) {

  auto print_map = [&os](const std::map<std::string, std::string> & properties) {
    unsigned int i=0;
    os << "[ ";
    for (const std::pair<std::string, std::string> & prop_and_val :  properties) {
      os << prop_and_val.first << "=\"" << prop_and_val.second << "\"";
      
      if (i != properties.size()-1)
	os << ", ";
      
      ++i;
    }
    os << " ];" << std::endl;
  };
  
  os << "graph {" << std::endl;

  const double small=0.025;

  // Nodes:
  for (MessagePasser<VARIABLE_KEY>*mp : ig.message_passers) {
    os << "\t\"" << mp << "\" ";

    std::map<std::string, std::string> node_properties;
    // default color:
    node_properties["color"] = "gray";
    node_properties["style"] = "filled";

    const HUGINMessagePasser<VARIABLE_KEY>* hmp = dynamic_cast<const HUGINMessagePasser<VARIABLE_KEY>* >(mp);
    if (hmp != NULL) {
      // HUGIN:
      node_properties["shape"] = "box";
      node_properties["color"] = "cyan";
      node_properties["style"] = "filled";
    
      const HUGINMessagePasser<VARIABLE_KEY>* he = dynamic_cast<const Hyperedge<VARIABLE_KEY>* >(hmp);
      if (he != NULL) {
	// Hyperedge:
	node_properties["style"] = "filled";
	node_properties["shape"] = "square";
	node_properties["color"] = "red";
	node_properties["width"] = to_string(small);
	node_properties["height"] = to_string(small);
      }
      else {
	// Not Hyperedge:
	
	if (hmp->joint_posterior().dimension() > 0) {
	  // HUGIN prior:

	  const std::vector<std::string> & vars = hmp->joint_posterior().ordered_variables();
	  for (unsigned int i=0; i<vars.size(); ++i) {
	    node_properties["label"] += to_string(vars[i]);
	    if (i != vars.size()-1)
	      node_properties["label"] += ",";
	  }
	}
      }
    }
    else {
      // Not HUGIN:
      
      const ConvolutionTreeMessagePasser<VARIABLE_KEY>* ctmp = dynamic_cast<const ConvolutionTreeMessagePasser<VARIABLE_KEY>* >(mp);
      if (ctmp != NULL) {
	// Convolution tree

	node_properties["color"] = "green";
	node_properties["shape"] = "triangle";
	node_properties["width"] = to_string(small);
	node_properties["height"] = to_string(small);
      }
      else {
	const ConstantMultiplierMessagePasser<VARIABLE_KEY>* cmmp = dynamic_cast<const ConstantMultiplierMessagePasser<VARIABLE_KEY>* >(mp);
	if (cmmp != NULL)
	  node_properties["color"] = "violet";
	node_properties["shape"] = "diamond";
	node_properties["width"] = to_string(small);
	node_properties["height"] = to_string(small);

	node_properties["label"] = "";
	const Vector<double> & scale = cmmp->scale();
	for (unsigned int i=0; i<scale.size(); ++i) {
	  node_properties["label"] += to_string(scale[i]);
	  if (i != scale.size()-1)
	    node_properties["label"] += ",";
	}
      }
    }
   
    if (node_properties.find("label") == node_properties.end())
      // Space is necessary for blank label in order to permit us to
      // height and width
      node_properties["label"] = " ";
    else
      node_properties["fontsize"] = "48";
    print_map(node_properties);
  }
  os << std::endl;

  std::set<Edge<VARIABLE_KEY>* > visited_edges;
  
  // Edges:
  for (MessagePasser<VARIABLE_KEY>*mp : ig.message_passers) {
    for (unsigned long k=0; k<mp->number_edges(); ++k) {
      Edge<VARIABLE_KEY>*edge = mp->get_edge_out(k);

      if (visited_edges.find(edge) != visited_edges.end())
	continue;
      
      visited_edges.insert(edge);
      visited_edges.insert(edge->get_opposite_edge_ptr());

      std::map<std::string, std::string> edge_properties;

      edge_properties["fontsize"] = "32";
      
      os << "\t\"" << edge->source << "\" -- \"" << edge->dest << "\"";
      for (unsigned int i=0; i<edge->variables_ptr->size(); ++i) {
	edge_properties["label"] += to_string( (*edge->variables_ptr)[i] );
	if (i != edge->variables_ptr->size()-1)
	  edge_properties["label"] += ",";
      }

      /*
      // Edge color:
      // Color edges red when they've passed:
      if (edge->dest->edge_received(edge->dest_edge_index))
	edge_properties["color"] = "red";
      // Color edges green when they are elligible to pass:
      else if (mp->ready_to_send_message_ab_initio(k))
	edge_properties["color"] = "green";
      */

      print_map(edge_properties);
    }
  }

  os << "}" << std::endl;
}

template <typename VARIABLE_KEY>
void write_graph_to_dot_file(const InferenceGraph<VARIABLE_KEY> & ig, const std::string & fname) {
  std::ofstream fout(fname);
  graph_to_dot(ig, fout);
  fout.close();
}

#endif
