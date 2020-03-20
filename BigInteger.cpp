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

    friend BigInteger operator - (const BigInteger&);

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

    enum SumType {
        PLUS,
        MINUS
    };

    enum DivideType {
        DIV,
        MOD
    };

    std::vector<int> bigInt;
    Sign sign;
    static const int baseLength = 3;
    static const int base = 1000;
    static const int naiveMulLenLimit = 128;

    void DeleteLeadingZeroes();
    void CorrectZeroSign();
    void ToNormalState();
    void FromString(const std::string&);
    BigInteger BuildHalf(const int&, const int&) const;
    BigInteger& BaseShift(const int&);

    static BigInteger& Sum(BigInteger&, const BigInteger&, SumType);
    static Sign DefineSumSign(Sign&, Sign&, CompareType, SumType);
    static BigInteger NaiveMultiply(BigInteger&, const BigInteger&);
    static BigInteger Karatsuba(BigInteger&, const BigInteger&);
    static BigInteger& Divide(BigInteger&, const BigInteger&, DivideType);
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

BigInteger& BigInteger::Sum(BigInteger& first, const BigInteger& second, BigInteger::SumType operand) {
    BigInteger result;

    Sign firstSign = first.sign;
    Sign secondSign = second.sign;
    result.sign = DefineSumSign(firstSign, secondSign, CompareByAbs(first,second), operand);

    int firstSize = first.bigInt.size();
    int secondSize = second.bigInt.size();
    int size = std::max(firstSize, secondSize) + 1;
    result.bigInt.resize(size, 0);
    first.bigInt.resize(size, 0);

    int baseOverFlow = 0;
    int firstTermCf;
    int secondTermCf;
    if (firstSign != secondSign) {
        if (CompareByAbs(first, second) == GREATER) {
            firstTermCf = 1;
            secondTermCf = -1;
        } else {
            firstTermCf = -1;
            secondTermCf = -1;
        }
    } else {
        firstTermCf = secondTermCf = 1;
    }

    for (int i = 0; i < size; ++i) {
        int firstTerm = (i < firstSize) ? first.bigInt[i] : 0;
        int secondTerm = (i < secondSize) ? second.bigInt[i] : 0;
        result.bigInt[i] = firstTerm * firstTermCf + secondTerm * secondTermCf + baseOverFlow;
        baseOverFlow = 0;
        if (result.bigInt[i] >= base) {
            baseOverFlow = 1;
            result.bigInt[i] %= base;
        } else if (result.bigInt[i] < 0) {
            baseOverFlow = -1;
            result.bigInt[i] += base;
        }
    }

    first = result;
    first.ToNormalState();
    return first;
}

BigInteger::Sign BigInteger::DefineSumSign(Sign& firstSign, Sign& secondSign, CompareType compare,
                                           BigInteger::SumType operand) {
    secondSign = operand == MINUS ? (secondSign == POSITIVE ? NEGATIVE : POSITIVE) : secondSign;
    if (firstSign == secondSign) {
        return firstSign;
    }
    return compare == GREATER ? (firstSign == POSITIVE ? POSITIVE : NEGATIVE) :
           (firstSign == POSITIVE ? NEGATIVE : POSITIVE);
}

BigInteger BigInteger::NaiveMultiply(BigInteger& first, const BigInteger& second) {
    BigInteger result;
    result.sign = first.sign == second.sign ? POSITIVE : NEGATIVE;

    int firstSize = first.bigInt.size();
    int secondSize = second.bigInt.size();
    int size = firstSize + secondSize;
    result.bigInt.resize(size + 1, 0);
    int baseOverFlow = 0;

    for (int i = 0; i < firstSize; ++i) {
        baseOverFlow = 0;
        for (int j = 0; j < secondSize; ++j) {
            int firstTerm = first.bigInt[i];
            int secondTerm = second.bigInt[j];
            result.bigInt[i+j] += firstTerm * secondTerm + baseOverFlow;
            if (result.bigInt[i+j] >= base) {
                baseOverFlow = result.bigInt[i+j] / base;
                result.bigInt[i+j] %= base;
            }
        }
        result.bigInt[i+secondSize] += baseOverFlow;
    }
    result.ToNormalState();
    return result;
}

BigInteger BigInteger::Karatsuba(BigInteger& first, const BigInteger& second) {
    int size = std::max(first.bigInt.size(), second.bigInt.size());
    if (first == BigInteger() || second == BigInteger()) {
        return BigInteger();
    }
    if (size <= 128) {
        return NaiveMultiply(first, second);
    }
    size += size % 2;
    int halfSize = size / 2 + size % 2;
    BigInteger tmpFirst = first;
    BigInteger tmpSecond = second;
    tmpFirst.bigInt.resize(size, 0);
    tmpSecond.bigInt.resize(size, 0);

    BigInteger firstL = tmpFirst.BuildHalf(halfSize, size);
    BigInteger firstR = tmpFirst.BuildHalf(0, halfSize);

    BigInteger secondL = tmpSecond.BuildHalf(halfSize, size);
    BigInteger secondR = tmpSecond.BuildHalf(0, halfSize);

    BigInteger firstProduct = Karatsuba(firstL, secondL);
    BigInteger secondProduct = Karatsuba(firstR, secondR);
    BigInteger thirdProduct = Karatsuba(firstL += firstR, secondL += secondR) - firstProduct - secondProduct;


    return firstProduct.BaseShift(size) + secondProduct + thirdProduct.BaseShift(size - halfSize);
}

BigInteger BigInteger::BuildHalf(const int& start, const int& end) const {
    BigInteger result;
    result.sign = sign;
    result.bigInt.clear();
    for (int i = start; i < end; ++i) {
        result.bigInt.push_back(bigInt[i]);
    }
    return result;
}

BigInteger& BigInteger::BaseShift(const int& power) {
    BigInteger temporaryNumber;
    temporaryNumber.bigInt.clear();
    temporaryNumber.sign = sign;
    int size = bigInt.size() + power;
    for (int i = 0; i < size; ++i) {
        temporaryNumber.bigInt.push_back(i < power ? 0 : bigInt[i-power]);
    }
    *this = temporaryNumber;
    return *this;
}

BigInteger& BigInteger::Divide(BigInteger& first, const BigInteger& second, DivideType type) {
    if (second == 0) {
        std::cerr << "Division by zero" << std::endl;
    }
    BigInteger result;
    result.sign = first.sign == second.sign ? POSITIVE : NEGATIVE;

    BigInteger divident = first;
    divident.sign = POSITIVE;

    BigInteger divisor = second;
    divisor.sign = POSITIVE;

    BigInteger baseShift;

    while (divident >= divisor) {
        baseShift = 1;
        while (divident >= divisor * base) {
            divisor *= base;
            baseShift *= base;
        }

        while (divident >= divisor) {
            divident -= divisor;
            result += baseShift;
        }

        divisor = second;
        divisor.sign = POSITIVE;
    }
    divident.ToNormalState();
    result.ToNormalState();
    first = type == DIV ? result : divident;
    return first;
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

BigInteger& BigInteger::operator=(const BigInteger& number) {
    this->bigInt = number.bigInt;
    this->sign = number.sign;
    return  *this;
}

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

BigInteger& operator++(BigInteger& number) {
    return number += 1;
}

BigInteger& operator--(BigInteger& number) {
    return number -= 1;
}

BigInteger operator++(BigInteger& number, int) {
    BigInteger oldNumber = number;
    number += 1;
    return oldNumber;
}

BigInteger operator--(BigInteger& number, int) {
    BigInteger oldNumber = number;
    number -= 1;
    return oldNumber;
}

BigInteger operator-(const BigInteger& number) {
    BigInteger result = number;
    result.sign = number.sign == BigInteger::POSITIVE ? BigInteger::NEGATIVE : BigInteger::POSITIVE;
    result.ToNormalState();
    return result;
}

BigInteger& operator+=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Sum(first, second, BigInteger::PLUS);
}

BigInteger& operator-=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Sum(first, second, BigInteger::MINUS);
}

BigInteger& operator*=(BigInteger& first, const BigInteger& second) {
    first = BigInteger::Karatsuba(first, second);
    return first;
}

BigInteger& operator/=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Divide(first, second, BigInteger::DIV);
}

BigInteger& operator%=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Divide(first, second, BigInteger::MOD);
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    return result += second;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second)  {
    BigInteger result = first;
    return result -= second;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    return result *= second;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    return result /= second;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    return result %= second;
}

BigInteger::operator bool() const {
    return (*this) != 0;
}

int main() {
    BigInteger a,b,c;
    std::cin >> a;
    std::cin >> b;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << (a / b) << std::endl;
    std::cout << (a % b) << std::endl;

    return 0;
}
