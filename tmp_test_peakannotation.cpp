// Quick test to verify PeakAnnotation methods work correctly
#include <OpenMS/METADATA/PeptideHit.h>
#include <iostream>
#include <vector>

using namespace OpenMS;

int main() {
    // Test basic construction and operators
    PeptideHit::PeakAnnotation a1, a2;
    
    a1.annotation = "y3";
    a1.charge = 2;
    a1.mz = 300.5;
    a1.intensity = 1000.0;
    
    a2.annotation = "b4";
    a2.charge = 1;
    a2.mz = 250.3;
    a2.intensity = 500.0;
    
    // Test operator==
    bool equal = (a1 == a2);  // Should be false
    bool equal_self = (a1 == a1);  // Should be true
    
    // Test operator<
    bool less = (a2 < a1);  // Should be true (smaller mz)
    
    // Test writePeakAnnotationsString_
    std::vector<PeptideHit::PeakAnnotation> annotations;
    annotations.push_back(a1);
    annotations.push_back(a2);
    
    String result;
    PeptideHit::PeakAnnotation::writePeakAnnotationsString_(result, annotations);
    
    std::cout << "Test completed successfully!" << std::endl;
    std::cout << "Equal: " << equal << ", Equal self: " << equal_self << std::endl;
    std::cout << "Less: " << less << std::endl;
    std::cout << "Annotation string: " << result.c_str() << std::endl;
    
    return 0;
}