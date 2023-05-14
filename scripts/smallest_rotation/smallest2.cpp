#include <iostream>

using namespace std;

int min_cyc(const string &s)
{
    int n = s.size();
    int res = 0;
    for (int l = 0; l < n; )
    {
        res = l;
        int r = l, p = l + 1;
        for (; r < n; ++r, ++p) /// If there is such string found, then its length wont exceed |s|
        {
            char c = (p < n) ? s[p] : s[p - n]; /// to avoid modulo
            if (s[r] > c) break;
            if (s[r] < c) r = l - 1;
        }
        l = max(r, l + p - r); /// just skip those (s2 = sx + sx + ... + sx + sy) cases
    }
    return res;
}

int main()
{
    string s;
    cin >> s;
    cout << min_cyc(s);
    return 0;
}
