import numpy as np

class Stochastic_function:
    '''
    Возвращает значение функции (метод get_value) для переданного значения времени
    param
    -----
    split_point_cnt: int
            количество членов в разложении
    func: Funtion
            объект функционала, включаеющего различные промежуточные вычисления
    period: float
                длина отрезка разложения
    variance_lst: lst
                список ранее вычисленных дисперсий
    c_coef_matrix: np.array, shape=(point_cnt, point_cnt)
                матрица рассчитанных коэфициентов c
    func_val_matrix: np.array, shape=(t_points.shape[0], self.split_point_cnt)
                значения функций в разложении
    '''
    def __init__(self, point_cnt, period, func) -> None:
        '''
        param
        -----
        point_cnt: int
                количество членов в разложении
        period: float
                длина отрезка разложения
        func: Funtion
                объект функционала, включаеющего различные промежуточные вычисления
        variance_lst: lst
                список ранее вычисленных дисперсий
        c_coef_matrix: np.array, shape=(point_cnt, point_cnt)
                матрица рассчитанных коэфициентов c
        '''
        self.split_point_cnt = point_cnt
        self.period = period
        self.func = func

        self.variance_lst = [0 for _ in range(self.split_point_cnt)]
        self.c_coef_matrix = np.zeros((self.split_point_cnt, self.split_point_cnt))

    def get_value(self, t_points):
        '''
        Возвращает значения случайной функции на переданном интервале
        param
        -----
        t_points: lst
                временные точки, в которых необходимо получить значения функции
        '''
        self.func.set_pendulums(self.split_point_cnt)
        self.func_val_matrix = np.zeros((t_points.shape[0], self.split_point_cnt))
        result = np.zeros((1,))
        for n in range(1, self.split_point_cnt):
            result += self.calculate_stochastic_value(n) * self.calculate_function(n, t_points)
        return result
    
    def calculate_c_coef(self, n, m) -> float:
        '''
        Возвращает значение коэффициента c_nm
        param
        -----
        n: int
        
        m: int
        '''
        if m == 1:
            if self.c_coef_matrix[n - 1, 0] == 0:
                self.c_coef_matrix[n - 1, 0] = - self.func.get_corr_value(n, 1) / self.variance_lst[0]
            return self.c_coef_matrix[n - 1, 0]
        temp_sum = 0
        for lambda_ind in range(1, m):
            if self.c_coef_matrix[n - 1, lambda_ind - 1] == 0:
                self.c_coef_matrix[n - 1, lambda_ind - 1] = self.calculate_c_coef(n, lambda_ind)
            if self.c_coef_matrix[m - 1, lambda_ind - 1] == 0:
                self.c_coef_matrix[m - 1, lambda_ind - 1] = self.calculate_c_coef(m - 1, lambda_ind)
            temp_sum += self.c_coef_matrix[n - 1, lambda_ind - 1] \
                        * self.c_coef_matrix[m - 1, lambda_ind - 1] * self.variance_lst[lambda_ind - 1]
        self.c_coef_matrix[n - 1, m - 1] = 1 / self.variance_lst[m - 1] * (temp_sum - self.func.get_corr_value(n, m))
        return self.c_coef_matrix[n - 1, m - 1]

    def calculate_variance(self, n: int) -> float:
        '''
        Возвращает значение дисперсии D_n
        param
        -----
        n: int
        '''
        if n == 1:
            return self.func.get_corr_value(1, 1)
        temp_sum = 0
        for i in range(1, n):
            temp_sum += np.abs(self.calculate_c_coef(n, i))**2 * self.variance_lst[i - 1]
        self.variance_lst[n - 1] = self.func.get_corr_value(n, n) - temp_sum
        return self.variance_lst[n - 1]

    def calculate_stochastic_value(self, n):
        '''
        Возвращает значение случайной величины V_n
        param
        -----
        n: int
        '''
        m = 0 # мат. ожидание
        variance = self.calculate_variance(n)
        return np.random.normal(m, np.sqrt(variance), 1)[0]
    
    def calculate_function(self, n, t_points):
        '''
        Вычисляет значения ф-ии в разложении
        param
        -----
        n: int
                номер функции в разложении
        t_points: lst
                временные точки, в которых необходимо получить значения функции
        '''
        if n == 1:
            self.func_val_matrix[0, :] = self.func.z_1(t_points) / self.variance_lst[0]
            return self.func_val_matrix[0, :]
        elif n == 2:
            self.func_val_matrix[1, :] = self.func.z_2n(t_points, 1) / self.variance_lst[1]
            return self.func_val_matrix[1, :]
        elif n == 3:
            return ((self.c_coef_matrix[2, 0] * self.func_val_matrix[0, :] + self.func.z_2n1(t_points, 1)) \
                    / self.variance_lst[3])
        elif n % 2 == 0:
            temp_sum = 0
            for m in range(1, n // 2):
                temp_sum += self.c_coef_matrix[n - 1, 2 * m - 1] * self.func_val_matrix[2 * m - 1, :]
            self.func_val_matrix[n - 1, :] = 1 / self.variance_lst[n - 1] \
                                            * (temp_sum + self.func.z_2n(t_points, n // 2))
            return self.func_val_matrix[n - 1, :]
        else:
            temp_sum = 0
            for m in range(1, n // 2 + 1):
                temp_sum += self.c_coef_matrix[n - 1, 2 * m - 2] * self.func_val_matrix[2 * m - 2, :]
            self.func_val_matrix[n - 1, :] = 1 / self.variance_lst[n - 1] \
                                            * (temp_sum + self.func.z_2n1(t_points, n // 2))
            return self.func_val_matrix[n - 1, :]