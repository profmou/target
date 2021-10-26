# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ReferenceEllipsoid:
    """参考椭球参数

    成员：
        __csid 参考坐标系名称，默认值'wgs84'
        __a  长半轴，单位：米
        __f  扁率
        __b  短半轴 b = (1-f)*a
        __e  偏心率 e = sqrt(a^2 - b^2)/a

    默认值：
    wgs84
        a = 6378137.0
        f = 1.0/298.257223563

    cgcs2000
        a = 6378137.0
        f = 1.0/298.257222101
    """

    def __init__(self, name='wgs84', a=None, f=None):
        """初始化参数后不可修改

        @param [in] name: 参考坐标系名称，'wgs84'(默认值), 'cgcs2000'
        @param [in] a: 长半轴，单位：米，float > 0.0
        @param [in] f: 扁率，0.0 < float < 1.0
        @return ReferenceEllipsoid
        """
        if not isinstance(name, str):
            print("name不是字符串，使用默认参数'wgs84'")
            nm = 'wgs84'
            aa = 6378137.0
            ff = 1.0/298.257223563
        elif name == 'wgs84':
            nm = name
            aa = 6378137.0
            ff = 1.0/298.257223563
        elif name == 'cgcs2000':
            nm = name
            aa = 6378137.0
            ff = 1.0/298.257222101
        else:
            nm = name

            arr_a = np.array(a)
            if arr_a.shape == () and isinstance(arr_a[()], np.number):
                # 数值标量
                aa = arr_a[()]
            elif arr_a.shape == (1,) and isinstance(arr_a[0], np.number):
                # 长度为1的数值向量
                aa = arr_a[0]
            else:
                print("a不是数值，使用默认参数'wgs84'的长半轴")
                aa = 6378137.0

            if aa < 0.0:
                print("a小于0，使用默认参数'wgs84'的长半轴")
                aa = 6378137.0

            arr_f = np.array(f)
            if arr_f.shape == () and isinstance(arr_f[()], np.number):
                # 数值标量
                ff = arr_f[()]
            elif arr_f.shape == (1,) and isinstance(arr_f[0], np.number):
                # 长度为1的数值向量
                ff = arr_f[0]
            else:
                print("f不是数值，使用默认参数'wgs84'的扁率")
                ff = 1.0/298.257223563

            if ff <= 0.0 or ff >= 1.0:
                print("f取值不在(0,1)之间，使用默认参数'wgs84'的扁率")
                ff = 1.0/298.257223563

        aa = np.array(aa, dtype=np.float64)
        ff = np.array(ff, dtype=np.float64)
        bb = (1-ff)*aa
        ee = (aa**2-bb**2)**0.5/aa

        self.__csid = nm
        self.__a = aa
        self.__f = ff
        self.__b = bb
        self.__e = ee

    def csid(self):
        """参考坐标系名称"""
        return self.__csid

    def a(self):
        """长半轴"""
        return self.__a

    def f(self):
        """扁平率"""
        return self.__f

    def b(self):
        """短半轴"""
        return self.__b

    def e(self):
        """偏心率"""
        return self.__e

def lla2ecef(lla, re=ReferenceEllipsoid()):
    """纬经高转直角坐标

    @param [in] lla 纬度(单位：度)、经度(单位：度)、高度(单位：米)
    @param [in] re 参考椭球参数
    @return xyz 单位：米
    """
    xyz = None
    if not isinstance(re, ReferenceEllipsoid):
        print("参考椭球参数不是由ReferenceEllipsoid生成，使用默认参数'wgs84'")
        rr = ReferenceEllipsoid()
        a = rr.a()
        e = rr.e()
    else:
        a = re.a()
        e = re.e()

    arr_lla = np.array(lla)
    if arr_lla.shape != (3,):
        print("纬经高参数个数不等于3")
    elif not isinstance(arr_lla[0], np.number):
        print("纬经高数据不是数值")
    else:
        lat, lon, alt = arr_lla.astype(np.float64)
        # 角度转为弧度
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)

        N = a/np.sqrt(1-(e*np.sin(lat))**2)
        x = (N + alt)*np.cos(lat)*np.cos(lon)
        y = (N + alt)*np.cos(lat)*np.sin(lon)
        z = (N*(1-e**2)+alt)*np.sin(lat)
        xyz = np.array([x, y, z])

    return xyz

def geodetic2enu(lla1, lla0, re=ReferenceEllipsoid()):
    """计算东北天坐标

    @param [in] lla1 待求坐标点，纬度(单位：度)、经度(单位：度)、高度(单位：米)
    @param [in] lla0 东北天坐标系原点，纬度(单位：度)、经度(单位：度)、高度(单位：米)
    @param [in] re 参考椭球参数
    @return enu 东北天坐标，单位：米
    """
    enu = None
    xyz1 = lla2ecef(lla1, re)
    xyz0 = lla2ecef(lla0, re)

    if (xyz1 is not None) and (xyz0 is not None):
        vv_xyz = np.array([xyz1 - xyz0]).T
        lat, lon, _ = np.array(lla0, dtype=np.float64)
        lat = np.deg2rad(lat)
        lon = np.deg2rad(lon)
        S = np.array([
            [-np.sin(lon), np.cos(lon), 0],
            [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
            [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]
            ])
        enu = np.dot(S, vv_xyz)
        enu = enu.T[0]

    return enu

def stat_enu(lla0):
    """统计"""
    def calc(data):
        enu = geodetic2enu(data, lla0)
        return pd.Series(data=enu, index=list('enu'))
    return calc

def cep(data, p):
    """圆概率误差"""
    if p > 100:
        p = 100
    sort_data = np.sort(data)    # sort只进行升序排序，降序需要手动逆序[::-1]
    idx = int(np.ceil(len(data)*p/100)) - 1
    if idx < 0:
        idx = 0
    return sort_data[idx]

def draw_target_chart(data, stat_value, title='打靶图', unit='m'):
    '''
    获取折线图
    :param data: pandas.DataFrame 两列数据：第一列，偏东值；第二列，偏北值
    :param stat_value double,统计值，如CEP95，drms2等
    :param title: str 图标题
    :param unit: str 单位
    :return:str 生成的图片的全文件名
    '''
    plt.rcParams['font.sans-serif'] = ['Simhei']
    plt.rcParams['axes.unicode_minus'] = False
    plt.figure(figsize=(5, 5), dpi=320)
    full_title = title + '(m)'
    plt.title(full_title)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_position(('data',0))
    ax.spines['bottom'].set_position(('data',0))
    xeast = data.iloc[:,[0]]
    ynorth = data.iloc[:,[1]]
    max_val = np.max([np.fabs(xeast.min()), np.fabs(xeast.max()), np.fabs(ynorth.min()), np.fabs(ynorth.max())])
    plt.xlim(-max_val - max_val*0.1, max_val + max_val*0.1)
    plt.ylim(-max_val - max_val*0.1, max_val + max_val*0.1)
    plt.text(0.52, 0.99, 'N', horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes)
    plt.text(0.99, 0.52, 'E', horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes)
    mid_val = max_val/2.0
    tickets = [-max_val, -mid_val, 0, mid_val, max_val]
    lbs = ['%.2f'%(-max_val), '%.2f'%(-mid_val), '0', '%.2f'%(mid_val), '%.2f'%(max_val)]
    plt.xticks(tickets, labels=lbs)
    del tickets[2]
    del lbs[2]
    plt.yticks(tickets)
    params = np.linspace(0, 2*np.pi, 100)
    max_x = np.cos(params)*max_val
    max_y = np.sin(params)*max_val
    plt.plot(max_x, max_y, color='black', linewidth=0.6)
    mid_x = np.cos(params)*mid_val
    mid_y = np.sin(params)*mid_val
    plt.plot(mid_x, mid_y, color='black', linewidth=0.6, linestyle='--')
    cep_x = np.cos(params)*stat_value
    cep_y = np.sin(params)*stat_value
    plt.plot(cep_x, cep_y, color='blue', linewidth=0.6)
    plt.annotate(s='CEP95:%.2f'%(stat_value), xy=(np.cos(np.pi/4.0)*stat_value,np.sin(np.pi/4.0)*stat_value),
                 xycoords='data', xytext=(mid_val, max_val),
                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=.2'))
    plt.scatter(xeast, ynorth, color='black', s=0.5)
    plt.show()

if __name__ == '__main__':
    df = pd.read_csv('lla.txt')
    lla0 = [df['lat'].mean(), df['lon'].mean(), df['alt'].mean()]
    temp = df.apply(stat_enu(lla0), axis=1)
    hers = (temp.e*temp.e + temp.n*temp.n)**0.5
    cep95 = cep(hers, 95)
    draw_target_chart(temp, cep95)
