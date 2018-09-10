#ifndef _SORT_H_
#define _SORT_H_

class Sort{
public:
    static void bubble(int* data, int n)
    {
        int i, j, temp;
        for (i = n - 1; i > 0; i--)
        {
            for (j = 0; j <= i - 1; j++)
            {
                if (data[j] > data[j + 1])
                {
                    temp = data[j];
                    data[j] = data[j + 1];
                    data[j + 1] = temp;
                }
            }
        }
    }
};

#endif
