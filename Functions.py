import numpy as np

class Functions:
    def __init__(self, period, *args) -> None:
        '''
        param
        -----
        period: float
                длина отрезка разложения
        args: lst
                параметры, использующиеся в корр. ф-ии (alpha, D)
        '''
        self.period = period
        self.alpha, self.D = args

    def z_1(self, t) -> float:
        '''
        Вычисляет значение промежуточной функции z_1(t)
        param
        -----
        t: float
                время
        '''
        result = []
        for i in range(t.shape[0]):
            result.append(self.D / self.alpha * (2 - np.exp(-self.alpha * t[i]) \
                - np.exp(-self.alpha*(self.period - t[i]))))
        return np.array(result)

    def z_2n(self, t, n) -> float:
        '''
        Вычисляет значение промежуточной функции z_2n(t)
        param
        -----
        t: float
                время
        n: int
                номер члена разложения
        '''
        result = []
        for i in range(t.shape[0]):
            result.append(self.D * self.alpha / (self.alpha**2 + self.freq_lst[n - 1]**2) \
                          *(2 * np.sin(self.freq_lst[n - 1] * t[i]) \
                            +self.freq_lst[n - 1]/self.alpha*np.exp(-self.alpha \
                            *t[i]) - np.exp(-self.alpha*(self.period-t[i])) \
                            *(np.sin(self.freq_lst[n - 1] * self.period) \
                            +self.freq_lst[n - 1]/self.alpha*np.cos(self.freq_lst[n - 1] * self.period))))
        return result
    
    def z_2n1(self, t, n) -> float:
        '''
        Вычисляет значение промежуточной функции z_2n+1(t)
        param
        -----
        t: float
                время
        n: int
                номер члена разложения
        '''
        result = []
        for i in range(t.shape[0]):
            result.append(self.D * self.alpha / (self.alpha**2 + self.freq_lst[n - 1]**2) \
                          *(2 * np.cos(self.freq_lst[n - 1] * t[i]) \
                            -np.exp(-self.alpha \
                            *t[i]) - np.exp(-self.alpha*(self.period-t[i])) \
                            *(np.cos(self.freq_lst[n - 1] * self.period) \
                            -self.freq_lst[n - 1]/self.alpha*np.sin(self.freq_lst[n - 1] * self.period))))

    def set_pendulums(self, point_cnt):
        '''
        Устанавливает частоты разложения
        param
        -----
        point_cnt: int
                количество членов в разложении
        '''
        self.freq_lst = [2 * np.pi * k / self.period for k in range(1, point_cnt)]

    def get_corr_value(self, n, m) -> float:
        '''
        Вычисляет значение коэффициента корреляции
        param
        -----
        n: int
        m: int
        '''
        if n == m == 1:
            return (2 * self.D / (self.alpha)**2 \
                    * (np.exp(-self.alpha * self.period) \
                    - 1 + self.alpha * self.period))
        elif (n == 1 and m % 2 == 0) or (n % 2 == 0 and m == 1):
            return 0
        elif (n == 1 and m % 2 != 0) or (n % 2 != 0 and m == 1):
            return ((2 * self.D * (1 - np.exp(-self.alpha*self.period))) \
                    / (self.alpha**2 + self.freq_lst[max(n, m) // 2]**2))
        elif n % 2 != 0 and m % 2 != 0 and n == m:
            return (self.D * self.alpha / (self.alpha**2 + self.freq_lst[max(n, m) // 2]**2) \
                    * (self.period - 2*self.alpha/(self.alpha**2 + self.freq_lst[max(n, m) // 2]**2) \
                    * (1 - np.exp(-self.alpha*self.period))))
        elif n % 2 != 0 and m % 2 != 0 and n != m:
            return (- 2 * self.D * self.alpha**2 * (1 - np.exp(-self.alpha * self.period)) \
                    / ((self.alpha**2 + self.freq_lst[n // 2]**2) \
                       *(self.alpha**2 + self.freq_lst[m // 2]**2)))
        elif n % 2 == 0 and m % 2 == 0 and n == m:
            return (self.D * self.alpha / (self.alpha**2 + self.freq_lst[max(n, m) // 2]**2) \
                    * (self.period + 2*self.freq_lst[max(n, m) // 2]**2 \
                       /(self.alpha*(self.alpha**2 + self.freq_lst[max(n, m) // 2]**2)) \
                    * (1 - np.exp(-self.alpha*self.period))))
        elif n % 2 == 0 and m % 2 == 0 and n != m:
            return (2 * self.D * self.freq_lst[n // 2] * self.freq_lst[m // 2] \
                    * (1 - np.exp(-self.alpha * self.period)) \
                    / ((self.alpha**2 + self.freq_lst[n // 2]**2) \
                       *(self.alpha**2 + self.freq_lst[m // 2]**2)))
        else:
            return 0
        