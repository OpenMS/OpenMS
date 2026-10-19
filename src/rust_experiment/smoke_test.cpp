#include <cstddef>
#include <cstring>

extern "C" std::ptrdiff_t openms_rust_count_elements(const unsigned char* data, std::size_t len);

int main()
{
  constexpr char xml[] = "<root><child/></root>";
  if (openms_rust_count_elements(reinterpret_cast<const unsigned char*>(xml), std::strlen(xml)) != 2)
  {
    return 1;
  }

  constexpr char malformed[] = "<root>";
  if (openms_rust_count_elements(reinterpret_cast<const unsigned char*>(malformed), std::strlen(malformed)) != -1)
  {
    return 1;
  }

  constexpr char mismatched[] = "<root><child></root>";
  return openms_rust_count_elements(reinterpret_cast<const unsigned char*>(mismatched), std::strlen(mismatched)) == -1 ? 0 : 1;
}
