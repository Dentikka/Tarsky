from sympy import *


def algorithm(pol):
    formula = ''
    max_deg = deg(pol)
    container = all_polynoms(pol)
    container = unique_polynoms(container)
    if pol.domain == ZZ:
        alg_with_const(pol, container, max_deg)
    else:
        container = drop_const(container)
        for i in range(3**len(container[0])):
            table = select_sign(container[0], i)
            for current_deg in range(1, max_deg + 1):
                for number_of_polynom in range(len(container[current_deg])):
                    current_polynom = container[current_deg][number_of_polynom]
                    add_row(table, current_polynom, container)
                    current_cell_number = 0
                    len_row = len(table[-1])
                    while current_cell_number < len_row - 1:
                        if table[-1][current_cell_number] * table[-1][current_cell_number + 1] == -1:
                            add_column(table, current_cell_number)
                            current_cell_number += 1
                            len_row += 1
                        current_cell_number += 1
            formula = ''
            formula = add_to_formula(table, container, formula)
            if formula == '':
                print('0 or')
            print(formula)
            table.clear()
        # formula = formula[:-4]
    container.clear()


def deg(pol):
    deg = degree(pol)
    if deg >= 0:
        return deg
    else:
        return 0


def norm(pol):
    return poly(str(pol.expr), symbols('x'))


def all_polynoms_prime(pol):
    container = [[] for _ in range(deg(pol) + 1)]
    container[deg(pol)].append(pol)

    for current_deg in range(deg(pol), 0, -1):
        for current_polynom in container[current_deg]:
            container[current_deg - 1].append(norm(diff(current_polynom)))
            remove_LT = current_polynom - LT(current_polynom)
            container[deg(remove_LT)].append(norm(remove_LT))
            container[0].append(poly(str(LC(current_polynom)), symbols('x')))

    return container


def all_polynoms_prime_c(pol, container):
    container[deg(pol)].append(pol)

    for current_deg in range(deg(pol), 0, -1):
        for current_polynom in container[current_deg]:
            container[current_deg - 1].append(norm(diff(current_polynom)))
            remove_LT = current_polynom - LT(current_polynom)
            container[deg(remove_LT)].append(norm(remove_LT))
            container[0].append(poly(str(LC(current_polynom)), symbols('x')))


def all_polynoms(pol):
    container = all_polynoms_prime(pol)

    for divisor_deg in range(deg(pol) - 1, 0, -1):
        for dividend_deg in range(deg(pol), divisor_deg - 1, -1):
            for divisor in container[divisor_deg]:
                for dividend in container[dividend_deg]:
                    all_polynoms_prime_c(norm(rem(dividend * LC(divisor) ** deg(dividend), divisor)), container)
    return container


def unique_polynoms(container):
    num_of_degs = len(container)
    for current_deg in range(num_of_degs):
        for current_polynom in container[current_deg]:
            count = container[current_deg].count(current_polynom)
            for i in range(count - 1):
                container[current_deg].remove(current_polynom)
    return container


def drop_const(container):
    cont = container[0]
    i = 0
    while i < len(cont):
        i += 1
        current_polynom = cont[i - 1]
        if ((current_polynom.expr < 0) == True) or ((current_polynom.expr > 0) == True) or ((current_polynom.expr == 0) == True):
            cont.remove(current_polynom)
            i -= 1
    return container


def to_trinary(x, number_of_digits):
    num = x
    base = 3
    newNum = ''
    for digit_number in range(number_of_digits):
        newNum = str(num % base) + newNum
        num //= base
    return newNum


def select_sign(constants, i):
    number_of_constants = len(constants)
    num_str = to_trinary(i, number_of_constants)
    table = [[0, 0], ]
    for const_num in range(len(constants)):
        const = constants[const_num]
        table.append([int(num_str[const_num]) - 1, int(num_str[const_num]) - 1])
    return table


def put_edge(table, current_polynom, container):
    deg_polynom = deg(current_polynom)
    oldest = poly(str(LC(current_polynom)), symbols('x'))
    if ((oldest.expr < 0) == True) or ((oldest.expr > 0) == True) or ((oldest.expr == 0) == True):
        if oldest.expr > 0:
            right_sign = 1
        elif oldest.expr < 0:
            right_sign = -1
        else:
            print('error')
    else:
        LC_index = container[0].index(oldest)
        right_sign = table[LC_index + 1][0]
        if right_sign == 0:
            remove_LT = norm(current_polynom - LT(current_polynom))
            if ((remove_LT.expr < 0) == True) or ((remove_LT.expr > 0) == True) or ((remove_LT.expr == 0) == True):
                if remove_LT.expr > 0:
                    right_sign = 1
                elif remove_LT.expr < 0:
                    right_sign = -1
                else:
                    right_sign = 0
            else:
                old_row_number = find_row(remove_LT, container)
                right_sign = table[old_row_number][-1]
            deg_polynom = deg(remove_LT)
    left_sign = right_sign
    if deg_polynom % 2:
        left_sign = -right_sign
    table[-1][0] = left_sign
    table[-1][-1] = right_sign


def find_polynom(row_number, container):
    counter = row_number - 1
    current_deg = 0
    while counter >= 0:
        counter -= len(container[current_deg])
        current_deg += 1
    current_deg -= 1
    counter += len(container[current_deg])
    return container[current_deg][counter]


def find_row(polynom, container):
    counter = 0
    polynom_degree = deg(polynom)
    for i in range(polynom_degree):
        counter += len(container[i])
    counter += container[polynom_degree].index(polynom)
    return counter + 1


def add_row(table, current_polynom, container):
    table_width = len(table[0])
    table_height = len(table)
    table.append([0 for _ in range(table_width)])
    table_height += 1
    put_edge(table, current_polynom, container)
    for column_number in range(1, table_width - 1):
        for row_number in range(1, table_height - 1):
            if table[row_number][column_number] == 0:
                divisor = find_polynom(row_number, container)
                dividend = current_polynom
                remainder = norm(rem(dividend * LC(divisor)**deg(dividend), divisor))
                if ((remainder.expr < 0) == True) or ((remainder.expr > 0) == True) or ((remainder.expr == 0) == True):
                    if remainder.expr > 0:
                        table[table_height - 1][column_number] = 1
                    elif remainder.expr < 0:
                        table[table_height - 1][column_number] = -1
                    else:
                        table[table_height - 1][column_number] = 0
                else:
                    remainder_row = find_row(remainder, container)
                    table[table_height - 1][column_number] = table[remainder_row][column_number]
    return table


def add_column(table, current_cell_number):
    table[-1].insert(current_cell_number + 1, 0)
    for current_row in table[:-1]:
        left_value = current_row[current_cell_number]
        right_value = current_row[current_cell_number + 1]
        if left_value != 0:
            current_row.insert(current_cell_number + 1, left_value)
        elif right_value != 0:
            current_row.insert(current_cell_number + 1, right_value)
        else:
            current_row.insert(current_cell_number + 1, 0)
    return table


def add_to_formula(table, container, formula):
    if (table[-1].count(0) > 0) and (len(container[0]) > 0):
        formula += '('
        for i in range(len(container[0])):
            current_const = container[0][i]
            if table[i + 1][0] == 1:
                formula += str(current_const.expr) + ' > 0 and '
            elif table[i + 1][0] == 0:
                formula += str(current_const.expr) + ' = 0 and '
            else:
                formula += str(current_const.expr) + ' < 0 and '
        formula = formula[:-5]
        formula += ') or '
    return formula


def sign(num):
    if num > 0:
        return 1
    elif num == 0:
        return 0
    else:
        return -1


def drop_zero(constants):
    count = constants.count(poly('0', symbols('x')))
    for _ in range(count):
        constants.remove(poly('0', symbols('x')))


def push_const(constants):
    number_of_constants = len(constants)
    table = [[0, 0], ]
    for const_num in range(len(constants)):
        const = constants[const_num].expr
        table.append([sign(const), sign(const)])
    return table


def alg_with_const(pol, container, max_deg):
    if pol != poly('0', symbols('x')):
        drop_zero(container[0])
    table = push_const(container[0])
    for current_deg in range(1, max_deg + 1):
        for number_of_polynom in range(len(container[current_deg])):
            current_polynom = container[current_deg][number_of_polynom]
            add_row(table, current_polynom, container)
            current_cell_number = 0
            len_row = len(table[-1])
            while current_cell_number < len_row - 1:
                if table[-1][current_cell_number] * table[-1][current_cell_number + 1] == -1:
                    add_column(table, current_cell_number)
                    current_cell_number += 1
                    len_row += 1
                current_cell_number += 1
    if (table[-1].count(0) > 0) and (len(container[0]) > 0):
        formula = '1'
    else:
        formula = '0'
    table.clear()
    print(formula)


if __name__ == '__main__':
    expr = str(input())
    pol = poly(expr, symbols('x'))
    algorithm(pol)