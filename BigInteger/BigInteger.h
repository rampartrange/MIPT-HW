#include <iostream>
#include <string>
#include <vector>

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

class BigInteger {
public:
    BigInteger();
    BigInteger(int);
    BigInteger(std::string);

    friend class Rational;

    friend std::ostream& operator<<(std::ostream&, const BigInteger&);
    friend std::istream& operator>>(std::istream&, BigInteger&);

    friend bool operator==(const BigInteger&, const BigInteger&);
    friend bool operator!=(const BigInteger&, const BigInteger&);
    friend bool operator>(const BigInteger&, const BigInteger&);
    friend bool operator<(const BigInteger&, const BigInteger&);
    friend bool operator>=(const BigInteger&, const BigInteger&);
    friend bool operator<=(const BigInteger&, const BigInteger&);

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

    std::string toString() const;

    friend BigInteger GCD(BigInteger, BigInteger);

private:
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
    std::string DecimalToString() const;

    static BigInteger& Sum(BigInteger&, const BigInteger&, SumType);
    static Sign DefineSumSign(Sign&, Sign&, CompareType, SumType);
    static BigInteger& NaiveMultiply(BigInteger&, const BigInteger&);
    static BigInteger& Karatsuba(BigInteger&, const BigInteger&);
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

BigInteger::BigInteger(std::string inStr) {
    FromString(inStr);
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

std::string BigInteger::toString() const {
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

std::string BigInteger::DecimalToString() const {
    std::string outStr;
    if (sign == NEGATIVE) {
        outStr += '-';
    }
    int size = bigInt.size();
    for (int i = size - 1; i >= 0; --i) {
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

BigInteger& BigInteger::Sum(BigInteger& first, const BigInteger& second, SumType operand) {
    BigInteger result;

    Sign firstSign = first.sign;
    Sign secondSign = second.sign;
    result.sign = DefineSumSign(firstSign, secondSign, CompareByAbs(first,second), operand);

    int firstSize = first.bigInt.size();
    int secondSize = second.bigInt.size();
    int size = std::max(firstSize, secondSize) + 1;
    result.bigInt.resize(size, 0);

    int baseOverFlow = 0;
    int firstTermCf;
    int secondTermCf;
    if (firstSign != secondSign) {
        if (CompareByAbs(first, second) == GREATER) {
            firstTermCf = 1;
            secondTermCf = -1;
        } else {
            firstTermCf = -1;
            secondTermCf = 1;
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

Sign BigInteger::DefineSumSign(Sign& firstSign, Sign& secondSign, CompareType compare,
                                           SumType operand) {
    secondSign = operand == MINUS ? (secondSign == POSITIVE ? NEGATIVE : POSITIVE) : secondSign;
    if (firstSign == secondSign) {
        return firstSign;
    }
    return compare == GREATER ? (firstSign == POSITIVE ? POSITIVE : NEGATIVE) :
           (firstSign == POSITIVE ? NEGATIVE : POSITIVE);
}

BigInteger& BigInteger::NaiveMultiply(BigInteger& first, const BigInteger& second) {
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
            baseOverFlow = result.bigInt[i+j] / base;
            result.bigInt[i+j] %= base;
        }
        result.bigInt[i+secondSize] += baseOverFlow;
    }
    result.ToNormalState();
    first = result;
    return first;
}

BigInteger& BigInteger::Karatsuba(BigInteger& first, const BigInteger& second) {
    int size = std::max(first.bigInt.size(), second.bigInt.size());

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

    BigInteger firstProduct = firstL;
    firstProduct = Karatsuba(firstProduct, secondL);
    BigInteger secondProduct = firstR;
    secondProduct = Karatsuba(secondProduct, secondR);
    BigInteger thirdProduct = Karatsuba(firstL += firstR, secondL += secondR) - firstProduct - secondProduct;

    first = (firstProduct.BaseShift(size) += secondProduct += thirdProduct.BaseShift(size-halfSize));
    //return firstProduct.BaseShift(size) + secondProduct + thirdProduct.BaseShift(size - halfSize);
    return first;
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
    Sign resultSign = first.sign == second.sign ? POSITIVE : NEGATIVE;
    Sign dividentSign = first.sign;

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
    result.sign = resultSign;
    divident.sign = dividentSign;
    divident.ToNormalState();
    result.ToNormalState();
    first = type == DIV ? result : divident;
    return first;
}

CompareType BigInteger::Compare(const BigInteger& first, const BigInteger& second) {
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

CompareType BigInteger::CompareByAbs(const BigInteger& first, const BigInteger& second) {
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
    out << number.toString();
    return out;

}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    std::string inStr;
    in >> inStr;
    number.FromString(inStr);
    return in;

}

bool operator==(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == EQUAL;
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) != EQUAL;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == GREATER;
}

bool operator<(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == LOWER;
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == EQUAL ||
           BigInteger::Compare(first, second) == GREATER;
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
    return BigInteger::Compare(first, second) == EQUAL ||
           BigInteger::Compare(first, second) == LOWER;
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
    result.sign = number.sign == POSITIVE ? NEGATIVE : POSITIVE;
    result.ToNormalState();
    return result;
}

BigInteger& operator+=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Sum(first, second, PLUS);
}

BigInteger& operator-=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Sum(first, second, MINUS);
}

BigInteger& operator*=(BigInteger& first, const BigInteger& second) {
    first = BigInteger::Karatsuba(first, second);
    return first;
}

BigInteger& operator/=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Divide(first, second, DIV);
}

BigInteger& operator%=(BigInteger& first, const BigInteger& second) {
    return BigInteger::Divide(first, second, MOD);
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

void Swap(BigInteger& first, BigInteger& second) {
    BigInteger swap = first;
    first = second;
    second = swap;
}

BigInteger GCD(BigInteger first, BigInteger second) {
    while (second != 0) {
        first %= second;
        Swap(first, second);
    }
    return first;
}

//RATIONAL

class Rational {
public:
    Rational();
    Rational(int);
    Rational(BigInteger);
    Rational(BigInteger, BigInteger);

    friend std::ostream& operator<<(std::ostream&, const Rational&);

    friend bool operator==(const Rational&, const Rational&);
    friend bool operator!=(const Rational&, const Rational&);
    friend bool operator>(const Rational&, const Rational&);
    friend bool operator<(const Rational&, const Rational&);
    friend bool operator>=(const Rational&, const Rational&);
    friend bool operator<=(const Rational&, const Rational&);

    Rational& operator=(const Rational&);

    friend Rational& operator++(Rational&);
    friend Rational& operator--(Rational&);
    friend Rational operator++(Rational&, int);
    friend Rational operator--(Rational&, int);

    friend Rational operator - (const Rational&);

    friend Rational& operator+=(Rational&, const Rational&);
    friend Rational& operator-=(Rational&, const Rational&);
    friend Rational& operator*=(Rational&, const Rational&);
    friend Rational& operator/=(Rational&, const Rational&);

    friend Rational operator+(const Rational&, const Rational&);
    friend Rational operator-(const Rational&, const Rational&);
    friend Rational operator*(const Rational&, const Rational&);
    friend Rational operator/(const Rational&, const Rational&);

    std::string toString() const;
    std::string asDecimal(size_t) const;

    explicit operator double();

private:


    BigInteger denominator;
    BigInteger numerator;
    Sign sign;

    void ToNormalState();
    void DecimalDivision(std::string&, BigInteger&, BigInteger&, size_t) const;

    static Rational& Sum(Rational&, const Rational&, SumType);
    static Sign DefineSumSign(Sign&, Sign&, CompareType, SumType);
    static Rational& Multiply(Rational&, const Rational&);
    static Rational& Divide(Rational&, const Rational&);
    static CompareType Compare(const Rational&, const Rational&);
    static CompareType CompareByAbs(const Rational&, const Rational&);

};

Rational::Rational() {
    denominator = 1;
    numerator = 0;
    sign = POSITIVE;
}

Rational::Rational(int number) {
    denominator = 1;
    numerator = number >= 0 ? number : -number;
    sign = number >= 0 ? POSITIVE : NEGATIVE;
}

Rational::Rational(BigInteger number) {
    denominator = 1;
    numerator = number >= 0 ? number : -number;
    sign = number.sign;
}

Rational::Rational(BigInteger num, BigInteger denom) {
    denominator = denom >= 0 ? denom : -denom;
    numerator = num >= 0 ? num : -num;
    sign = num.sign == denom.sign ? POSITIVE : NEGATIVE;
}


std::string Rational::toString() const {
    return (sign == POSITIVE ? "" : "-") + numerator.toString() +
           (denominator == 1 ? "" : ("/" + denominator.toString()));
}

std::string Rational::asDecimal(size_t precision = 0) const{

    std::string output = sign == POSITIVE ? "" : "-";

    output += (numerator / denominator).toString();

    BigInteger num = numerator;
    BigInteger den = denominator;
    if (precision > 0) {
        output += '.';
        std::string decimalOutput = "";
        for (size_t i = 0; i < precision; ++i){
            num %= den;
            num *= BigInteger::base;
            decimalOutput += (num / den).DecimalToString();
        }
        for (size_t i = 0; i < precision; ++i) {
            output += decimalOutput[i];
        }
    }
    return output;
}

void Rational::ToNormalState() {
    if (numerator == 0) {
        sign = POSITIVE;
        denominator = 1;
        return;
    }

    BigInteger gcd = GCD(numerator, denominator);
    numerator /= gcd;
    denominator /= gcd;
    numerator.sign = POSITIVE;
    denominator.sign = POSITIVE;
}

void Rational::DecimalDivision(std::string& output, BigInteger& num, BigInteger& den,  size_t precision) const{
    BigInteger baseShift;
    BigInteger result;
    for (size_t i = 0; i < precision; ++i) {
        baseShift = 1;

        result = 0;
        num *= BigInteger::base;

        while (num >= den * BigInteger::base) {
            den *= BigInteger::base;
            baseShift *= BigInteger::base;
        }

        while (num >= den) {
            num -= den;
            result += baseShift;
        }

        output += result.toString();


        den = denominator;
    }
}

Rational& Rational::Sum(Rational& first, const Rational& second, SumType operand) {
    Rational result;
    Sign firstSign = first.sign;
    Sign secondSign = second.sign;
    result.sign = DefineSumSign(firstSign, secondSign, CompareByAbs(first, second), operand);

    int firstTermCf;
    int secondTermCf;
    if (firstSign != secondSign) {
        if (CompareByAbs(first, second) == GREATER) {
            firstTermCf = 1;
            secondTermCf = -1;
        } else {
            firstTermCf = -1;
            secondTermCf = 1;
        }
    } else {
        firstTermCf = secondTermCf = 1;
    }

    result.numerator = first.numerator * second.denominator * firstTermCf;
    result.numerator += second.numerator * first.denominator * secondTermCf;
    result.denominator = first.denominator * second.denominator;
    result.ToNormalState();
    first = result;
    return first;

}

Sign Rational::DefineSumSign(Sign& firstSign, Sign& secondSign, CompareType compare,
                               SumType operand) {
    secondSign = operand == MINUS ? (secondSign == POSITIVE ? NEGATIVE : POSITIVE) : secondSign;
    if (firstSign == secondSign) {
        return firstSign;
    }
    return compare == GREATER ? (firstSign == POSITIVE ? POSITIVE : NEGATIVE) :
           (firstSign == POSITIVE ? NEGATIVE : POSITIVE);
}

Rational& Rational::Multiply(Rational& first, const Rational& second) {
    first.numerator *= second.numerator;
    first.denominator *= second.denominator;
    first.sign = first.sign == second.sign ? POSITIVE : NEGATIVE;
    first.ToNormalState();
    return first;
}

Rational& Rational::Divide(Rational& first, const Rational& second) {
    if (second == 0) {
        std::cerr << "Division by zero";
    }
    first.numerator *= second.denominator;
    first.denominator *= second.numerator;
    first.sign = first.sign == second.sign ? POSITIVE : NEGATIVE;
    first.ToNormalState();
    return first;
}

CompareType Rational::Compare(const Rational& first, const Rational& second) {
    if (first.sign == POSITIVE && second.sign == NEGATIVE) {
        return GREATER;
    }
    if (first.sign == NEGATIVE && second.sign == POSITIVE) {
        return LOWER;
    }

    return  (first.sign == POSITIVE ? CompareByAbs(first, second) : CompareByAbs(second, first));
}

CompareType Rational::CompareByAbs(const Rational& first, const Rational& second) {
    return BigInteger::CompareByAbs(first.numerator * second.denominator, second.numerator * first.denominator);
}

//Operator overload

std::ostream& operator<<(std::ostream& out, const Rational& number) {
    out << number.toString();
    return out;

}

bool operator==(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) == EQUAL;
}

bool operator!=(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) != EQUAL;
}

bool operator>(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) == GREATER;
}

bool operator<(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) == LOWER;
}

bool operator>=(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) == EQUAL ||
           Rational::Compare(first, second) == GREATER;
}

bool operator<=(const Rational& first, const Rational& second) {
    return Rational::Compare(first, second) == EQUAL ||
           Rational::Compare(first, second) == LOWER;
}


Rational& Rational::operator=(const Rational& number) {
    numerator = number.numerator;
    denominator = number.denominator;
    sign = number.sign;
    return  *this;
}

Rational& operator++(Rational& number) {
    return number += 1;
}

Rational& operator--(Rational& number) {
    return number -= 1;
}

Rational operator++(Rational& number, int) {
    Rational result = number;
    number += 1;
    return result;
}

Rational operator--(Rational& number, int) {
    Rational result = number;
    number -= 1;
    return result;
}

Rational operator-(const Rational& number) {
    Rational result = number;
    result.sign = number.sign == POSITIVE ? NEGATIVE : POSITIVE;
    return result;
}

Rational& operator+=(Rational& first, const Rational& second) {
    return Rational::Sum(first, second, PLUS);
}

Rational& operator-=(Rational& first, const Rational& second) {
    return Rational::Sum(first, second, MINUS);
}

Rational& operator*=(Rational& first, const Rational& second) {
    return Rational::Multiply(first, second);
}

Rational& operator/=(Rational& first, const Rational& second) {
    return Rational::Divide(first, second);
}

Rational operator+(const Rational& first, const Rational& second) {
    Rational result = first;
    return result += second;
}
Rational operator-(const Rational& first, const Rational& second) {
    Rational result = first;
    return result -= second;
}
Rational operator*(const Rational& first, const Rational& second) {
    Rational result = first;
    return result *= second;
}
Rational operator/(const Rational& first, const Rational& second) {
    Rational result = first;
    return result /= second;
}

Rational::operator double() {
    return std::stod(asDecimal(1024));
}
