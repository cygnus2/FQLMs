#include <iostream>
#include <vector>

int main()
{
    std::vector<char> str;
    str.reserve(6);
    for(int i=0;i<6;i++) str[i]='0';
      str[6]='\0';

    for (auto ii : str) {
            std::cout << ii << ' ';
    }
    std::cout << '\n';
}
