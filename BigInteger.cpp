#include <iostream>
#include <string>
#include <vector>

class BigInteger {
public:
    BigInteger();
    BigInteger(int);

    friend std::ostream& operator<<(std::ostream&, const BigInteger&);
    friend std::istream& operator>>(std::istream&, BigInteger&);

    friend bool operator==(const BigInteger &, const BigInteger &);
    friend bool operator!=(const BigInteger &, const BigInteger &);
    friend bool operator>(const BigInteger &, const BigInteger &);
    friend bool operator<(const BigInteger &, const BigInteger &);
    friend bool operator>=(const BigInteger &, const BigInteger &);
    friend bool operator<=(const BigInteger &, const BigInteger &);

    BigInteger& operator=(const BigInteger&);

    friend BigInteger& operator++(BigInteger&);
    friend BigInteger& operator--(BigInteger&);
    friend BigInteger operator++(BigInteger&, int);
    friend BigInteger operator--(BigInteger&, int);

    friend BigInteger& operator+=(BigInteger&, const BigInteger&);
    friend BigInteger& operator-=(BigInteger&, const BigInteger&);
    friend BigInteger& operator*=(BigInteger&, const BigInteger&);
    friend BigInteger& operator/=(BigInteger&, const BigInteger&);
    friend BigInteger& operator%=(BigInteger&, const BigInteger&);

    friend BigInteger operator+(const BigInteger&, const BigInteger&);
    friend BigInteger operator-(const BigInteger&, const BigInteger&);
    friend BigInteger operator*(const BigInteger&, const BigInteger&);
    friend BigInteger operator/(const BigInteger&, const BigInteger&);
    friend BigInteger operator%(const BigInteger&, const BigInteger&);

    explicit operator bool() const;

    std::string ToString() const;
    void CHECK(const BigInteger&, const BigInteger&);

private:
    enum Sign {
        POSITIVE = 0,
        NEGATIVE = 1
    };

    enum CompareType {
        LOWER,
        EQUAL,
        GREATER
    };

    enum OperandType {
        PLUS,
        MINUS
    };

    std::vector<int> bigInt;
    Sign sign;
    static const int baseLength = 3;
    static const int base = 1000;

    void DeleteLeadingZeroes();
    void CorrectZeroSign();
    void ToNormalState();
    void FromString(const std::string&);

    static BigInteger& Sum(BigInteger&, const BigInteger&, OperandType);
    static Sign DefineSumSign(BigInteger&, const BigInteger&, OperandType);
    static CompareType Compare(const BigInteger&, const BigInteger&);
    static CompareType CompareByAbs(const BigInteger&, const BigInteger&);
};

BigInteger::BigInteger() {
    bigInt.clear();
    bigInt.push_back(0);
    sign = POSITIVE;
}

BigInteger::BigInteger(int number) {
    bigInt.clear();
    sign = number >= 0 ? POSITIVE : NEGATIVE;
    number = abs(number);

    if (number == 0) {
        bigInt.push_back(0);
    }

    while (number > 0) {
        bigInt.push_back(number % base);
        number /= base;
    }
}

void BigInteger::DeleteLeadingZeroes() {
    while (bigInt.size() > 1 && bigInt.back() == 0) {
        bigInt.pop_back();
    }
}

void BigInteger::CorrectZeroSign() {
    if (bigInt.size() == 1 && bigInt[0] == 0) {
        sign = POSITIVE;
    }
}

void BigInteger::ToNormalState() {
    DeleteLeadingZeroes();
    CorrectZeroSign();
}

void BigInteger::FromString(const std::string& inStr) {
    bigInt.clear();

    if (inStr.length() == 0) {
        return;
    }
    sign = inStr[0] == '-' ? NEGATIVE : POSITIVE;
    int strSize = inStr.length();
    for (int i = strSize; i > int(sign); i -= baseLength) {
        if (i <= baseLength) {
            bigInt.push_back(std::stoi(inStr.substr(int(sign), i - int(sign))));
        } else {
            bigInt.push_back(std::stoi(inStr.substr(i - baseLength, baseLength)));
        }
    }
    ToNormalState();
}

std::string BigInteger::ToString() const {
    std::string outStr;
    if (sign == NEGATIVE) {
        outStr += '-';
    }
    int size = bigInt.size();
    outStr += std::to_string(bigInt[size - 1]);
    for (int i = size - 2; i >= 0; --i) {
        int additionalZeroes = 0;
        int leftB = 1;
        int rightB = 10;
        while (not(bigInt[i] >= base / rightB && bigInt[i] < base / leftB) && base / rightB > 1) {
            ++additionalZeroes;
            leftB *= 10;
            rightB *= 10;
        }
        for (int j = 0; j < additionalZeroes; ++j) {
            outStr += '0';
        }
        outStr += std::to_string(bigInt[i]);
    }
    return outStr;

}

BigInteger::CompareType BigInteger::Compare(const BigInteger& first, const BigInteger& second) {
    if (first.sign == POSITIVE && second.sign == NEGATIVE) {
        return GREATER;
    } else if (first.sign == NEGATIVE && second.sign == POSITIVE) {
        return LOWER;
    }
    if (first.sign == POSITIVE || CompareByAbs(first, second) == EQUAL) {
        return CompareByAbs(first, second);
    } else {
        return CompareByAbs(first, second) == GREATER ? LOWER : GREATER;
    }

}

BigInteger::CompareType BigInteger::CompareByAbs(const BigInteger& first, const BigInteger& second) {
    if (first.bigInt.size() > second.bigInt.size()) {
        return GREATER;
    } else if (first.bigInt.size() < second.bigInt.size()) {
        return LOWER;
    }
    int size = first.bigInt.size();
    for (int i = size - 1; i >= 0; --i) {
        if (first.bigInt[i] > second.bigInt[i]) {
            return GREATER;
        } else if (first.bigInt[i] < second.bigInt[i]) {
            return LOWER;
        }
    }
    return EQUAL;
}

//Operator overload
std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.ToString();
    return out;

}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    std::string inStr;
    in >> inStr;
    number.FromString(inStr);
    return in;

}

bool operator==(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == BigInteger::EQUAL;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) != BigInteger::EQUAL;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == BigInteger::GREATER;
}

bool operator<(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == BigInteger::LOWER;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == BigInteger::EQUAL ||
           BigInteger::Compare(first, second) == BigInteger::GREATER;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == BigInteger::EQUAL ||
           BigInteger::Compare(first, second) == BigInteger::LOWER;
}

void BigInteger::CHECK(const BigInteger& a, const BigInteger& b) {
    std::cout<<Compare(a, b);
}

int main() {
    BigInteger a,b;
    std::cin >> a;
    std::cin >> b;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << (a == b) << std::endl;
    std::cout << (a != b) << std::endl;
    std::cout << (a > b) << std::endl;
    std::cout << (a < b) << std::endl;
    std::cout << (a >= b) << std::endl;
    std::cout << (a <= b) << std::endl;
    return 0;
}
