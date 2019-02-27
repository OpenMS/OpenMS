#ifndef _Elements_HPP
#define _Elements_HPP

struct Isotope {
  std::string name;
  double mass;
  double abundance;
};

// To enable std::set<Isotope>
bool operator <(const Isotope & lhs, const Isotope & rhs) {
  return lhs.name < rhs.name || (lhs.name == rhs.name && lhs.mass < rhs.mass);
}

std::ostream & operator<<(std::ostream & os, const Isotope & rhs) {
  os << rhs.name << ": mass=" << rhs.mass << " abundance=" << rhs.abundance;
  return os;
}

class Elements {
protected:
  std::map<std::string, std::vector<Isotope> > _isotope_list;
  
public:
  
  Elements(const std::string & isotop_file) {
    std::ifstream myfile(isotop_file);
    assert(myfile.is_open() == true && "Error: File not found");
    
    std::string line;
    std::string element;
    double mass;
    double min_abundance;
    double max_abundance;
    
    while ( std::getline(myfile,line) ) {
      std::istringstream ist(line);
      ist >> element;
      ist >> mass;
      ist >> min_abundance;
      ist >> max_abundance;

      Isotope iso = {element, mass, (max_abundance + min_abundance)/2};
      _isotope_list[element].push_back(iso);
    }
    myfile.close();
  }
  
  void print_elements_list() const {
    std::cout << "[ ";
    for(auto const & key: _isotope_list) {
      std::cout << "[";
      
      for(unsigned long i=0; i+1<key.second.size(); ++i) {
        std::cout << key.second[i] << ", ";
      }
      std::cout << key.second.back() << "] ";
    }
    std::cout << "]" << std::endl;
  }
  
  std::map<std::string, std::vector<Isotope> >::const_iterator find(const std::string & key) const {
    return _isotope_list.find(key);
  }
  
  std::map<std::string, std::vector<Isotope> >::const_iterator begin() const {
    return _isotope_list.begin();
  }
  
  std::map<std::string, std::vector<Isotope> >::const_iterator end() const {
    return _isotope_list.end();
  }
  
  unsigned long size() const {
    return _isotope_list.size();
  }
  
  const std::vector<Isotope> get(const std::string & key) const {
    return _isotope_list.at(key);
  }
  
};

#endif
