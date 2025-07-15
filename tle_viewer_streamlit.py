import streamlit as st
import requests
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from io import StringIO

def validate_tle(name, line1, line2):
    """TLEデータの基本的な検証"""
    if not name or not line1 or not line2:
        return False
    
    # TLE行の長さチェック
    if len(line1) != 69 or len(line2) != 69:
        return False
    
    # 行番号チェック
    if not line1.startswith('1 ') or not line2.startswith('2 '):
        return False
    
    return True

def fetch_tle_data(url="https://celestrak.org/NORAD/elements/gp.php?GROUP=last-30-days&FORMAT=tle"):
    """CelestrakからTLEデータを取得"""
    try:
        st.info("TLEデータを取得中...")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        lines = response.text.strip().split('\n')
        tle_data = []
        
        for i in range(0, len(lines), 3):
            if i + 2 < len(lines):
                name = lines[i].strip()
                line1 = lines[i + 1].strip()
                line2 = lines[i + 2].strip()
                
                if validate_tle(name, line1, line2):
                    tle_data.append({
                        'name': name,
                        'line1': line1,
                        'line2': line2
                    })
        
        st.success(f"{len(tle_data)}個の衛星データを取得しました")
        return tle_data
        
    except Exception as e:
        st.error(f"TLEデータの取得に失敗しました: {e}")
        return None

def parse_tle(tle_entry):
    """TLEデータを解析して軌道要素を抽出"""
    line1 = tle_entry['line1']
    line2 = tle_entry['line2']
    
    try:
        # 軌道要素の解析
        epoch_year = int(line1[18:20])
        epoch_day = float(line1[20:32])
        inclination = float(line2[8:16])
        raan = float(line2[17:25])
        eccentricity = float('0.' + line2[26:33])
        arg_perigee = float(line2[34:42])
        mean_anomaly = float(line2[43:51])
        mean_motion = float(line2[52:63])
        
        # 軌道周期の計算
        period = 24 * 60 * 60 / mean_motion  # 秒単位
        
        return {
            'name': tle_entry['name'],
            'epoch_year': epoch_year,
            'epoch_day': epoch_day,
            'inclination': inclination,
            'raan': raan,
            'eccentricity': eccentricity,
            'arg_perigee': arg_perigee,
            'mean_anomaly': mean_anomaly,
            'mean_motion': mean_motion,
            'period': period
        }
        
    except Exception as e:
        st.error(f"TLE解析エラー ({tle_entry['name']}): {e}")
        return None

def calculate_orbit(orbit_elements, end_time=500):
    """3次元軌道計算"""
    # 定数の設定
    G = 6.67430e-11             # 万有引力定数 [m^3 kg^-1 s^-2]
    M_earth = 5.972e24          # 地球質量 [kg]
    re = 6371012.0              # 地球半径 [m]
    g = G * M_earth / (re * re) # 重力加速度 [m/s^2]
    dt = 120.0                  # 時刻計算間隔 [s]
    eps = 0.00001               # 収束精度
    
    # 軌道要素の取得
    n = orbit_elements['mean_motion'] * 2.0 * np.pi / (24.0 * 60.0 * 60.0)  # 平均運動 [rad/s]
    e0 = orbit_elements['eccentricity']                                     # 離心率
    i = np.radians(orbit_elements['inclination'])                           # 軌道傾斜角 [rad]
    Ω = np.radians(orbit_elements['raan'])                                  # 昇交点赤経 [rad]
    ω = np.radians(orbit_elements['arg_perigee'])                           # 近地点引数 [rad]
    M0 = np.radians(orbit_elements['mean_anomaly'])                         # 平均近点角 [rad]
    
    # 軌道長半径の計算
    a = ((g * re * re) / (n * n)) ** (1.0/3.0)  # [m]
    
    # 座標変換行列の作成
    R_ω = np.array([[np.cos(ω), -np.sin(ω), 0],
                    [np.sin(ω), np.cos(ω), 0],
                    [0, 0, 1]])
    
    R_i = np.array([[1, 0, 0],
                    [0, np.cos(i), -np.sin(i)],
                    [0, np.sin(i), np.cos(i)]])
    
    R_Ω = np.array([[np.cos(Ω), -np.sin(Ω), 0],
                    [np.sin(Ω), np.cos(Ω), 0],
                    [0, 0, 1]])
    
    positions = []
    velocities = []
    t = 0
    
    # 軌道計算
    for j in range(end_time):
        t = j * dt
        M = M0 + n * t
        
        Eds = M
        while True:
            Edn = Eds - ((Eds - e0 * np.sin(Eds) - M) / (1 - e0 * np.cos(Eds)))
            if abs(Edn - Eds) < eps:
                break
            Eds = Edn
            
        # 位置計算
        xp = (a * np.cos(Edn) - a * e0) / 1000.0  # [km]
        yp = a * np.sqrt(1 - e0*e0) * np.sin(Edn) / 1000.0  # [km]
        zp = 0.0  # [km]
        
        # 速度計算（軌道面上）
        vxp = -n * a * np.sin(Edn) / (1 - e0 * np.cos(Edn)) / 1000.0  # [km/s]
        vyp = n * a * np.sqrt(1 - e0*e0) * np.cos(Edn) / (1 - e0 * np.cos(Edn)) / 1000.0  # [km/s]
        vzp = 0.0  # [km/s]
        
        # 座標変換
        pos = np.array([xp, yp, zp])
        pos = R_Ω @ R_i @ R_ω @ pos
        
        vel = np.array([vxp, vyp, vzp])
        vel = R_Ω @ R_i @ R_ω @ vel
        
        positions.append(pos)
        velocities.append(vel)
        
    positions = np.array(positions).T
    velocities = np.array(velocities).T
    return positions[0], positions[1], positions[2], velocities[0], velocities[1], velocities[2]

def parse_tle_text(tle_text):
    """テキストからTLEデータを解析"""
    lines = tle_text.strip().split('\n')
    tle_data = []
    
    for i in range(0, len(lines), 3):
        if i + 2 < len(lines):
            name = lines[i].strip()
            line1 = lines[i + 1].strip()
            line2 = lines[i + 2].strip()
            
            if validate_tle(name, line1, line2):
                tle_data.append({
                    'name': name,
                    'line1': line1,
                    'line2': line2
                })
    
    return tle_data

def load_tle_from_file(file_path="tle_cache.txt"):
    """ローカルファイルからTLEデータを読み込む"""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            tle_text = f.read()
        tle_data = parse_tle_text(tle_text)
        return tle_data
    except Exception as e:
        st.error(f"ローカルTLEファイルの読み込みに失敗しました: {e}")
        return None

def create_orbit_plot(tle_data, end_time, sat_num=5):
    """軌道プロットを作成"""
    earth_radius = 6378.0  # 地球半径 [km]
    
    # 地球の描画
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x_earth = earth_radius * np.outer(np.cos(u), np.sin(v))
    y_earth = earth_radius * np.outer(np.sin(u), np.sin(v))
    z_earth = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))

    fig = go.Figure()

    # 地球
    fig.add_surface(x=x_earth, y=y_earth, z=z_earth, 
                    colorscale=[[0, 'lightblue'], [1, 'lightblue']],
                    opacity=0.3,
                    showscale=False,
                    surfacecolor=np.ones((100, 100)),
                    showlegend=False,
                    hoverinfo='skip')
    
    # 東京の座標: 緯度 35.6762°N, 経度 139.6503°E
    lat_tokyo = np.radians(35.6762)  # 緯度（ラジアン）
    lon_tokyo = np.radians(139.6503)  # 経度（ラジアン）

    # 地心座標系での東京の位置
    x_tokyo = earth_radius * np.cos(lat_tokyo) * np.cos(lon_tokyo)
    y_tokyo = earth_radius * np.cos(lat_tokyo) * np.sin(lon_tokyo)
    z_tokyo = earth_radius * np.sin(lat_tokyo)
    
    fig.add_scatter3d(x=[x_tokyo], y=[y_tokyo], z=[z_tokyo],
                      mode='markers',
                      name='東京',
                      marker=dict(size=1, color='red', symbol='circle'),
                      hovertemplate='東京<br>緯度: 35.6762°N<br>経度: 139.6503°E<extra></extra>')

    colors = ['red', 'green', 'blue', 'yellow', 'purple', 'orange', 'pink', 'brown', 'gray', 'black'] * (sat_num // 10 + 1)
    
    for i, tle_entry in enumerate(tle_data[:sat_num]):
        orbit_elements = parse_tle(tle_entry)
        if orbit_elements:
            x, y, z, vx, vy, vz = calculate_orbit(orbit_elements, end_time)
            
            hover_texts = []
            for j in range(len(x)):
                speed = np.sqrt(vx[j]**2 + vy[j]**2 + vz[j]**2)
                # 高度計算（地球中心からの距離 - 地球半径）
                altitude = np.sqrt(x[j]**2 + y[j]**2 + z[j]**2) - 6378.0  # km
                hover_text = f"衛星: {orbit_elements['name']}<br>"
                hover_text += f"位置: ({x[j]:.1f}, {y[j]:.1f}, {z[j]:.1f}) km<br>"
                hover_text += f"高度: {altitude:.1f} km<br>"
                hover_text += f"速度: {speed:.4f} km/s<br>"
                hover_text += f"速度成分: ({vx[j]:.2f}, {vy[j]:.2f}, {vz[j]:.2f}) km/s"
                hover_texts.append(hover_text)
            
            fig.add_scatter3d(x=x, y=y, z=z,
                            mode='lines',
                            name=orbit_elements['name'][:20],
                            line=dict(color=colors[i], width=2),
                            hovertemplate='%{text}<extra></extra>',
                            text=hover_texts)

    fig.update_layout(
        scene=dict(
            xaxis_title='X [km]',
            yaxis_title='Y [km]',
            zaxis_title='Z [km]',
            aspectmode='data'
        ),
        showlegend=True,
        height=1000
    )
    
    return fig

# --- 追加: 軌道長半径計算関数 ---
def calculate_semi_major_axis(mean_motion):
    # mean_motion: [rev/day] -> [rad/s]に変換して軌道長半径aを計算
    # 公式: a = (mu / (n^2))^(1/3)
    # mu = GM = 3.986004418e14 [m^3/s^2]
    # n = mean_motion * 2pi / 86400 [rad/s]
    mu = 3.986004418e14
    n = mean_motion * 2 * np.pi / 86400
    a = (mu / (n ** 2)) ** (1/3) / 1000  # [km]
    return a

def main():
    st.set_page_config(page_title="TLE Viewer", layout="wide")
    
    st.header("🚀 TLE Viewer - 衛星軌道可視化ツール")
    
    # サイドバーにコントロールを配置
    with st.sidebar:
        st.header("設定")
        
        # データソース選択
        data_source = st.radio(
            "データソース",
            ["TLEファイル", "手動入力"],
            help="TLEデータの取得方法を選択してください"
        )
        
        # オンライン取得の場合
        if data_source == "オンライン取得":
            url_options = {
                "過去30日間": "https://celestrak.org/NORAD/elements/gp.php?GROUP=last-30-days&FORMAT=tle",
                "気象衛星": "https://celestrak.org/NORAD/elements/gp.php?GROUP=weather&FORMAT=tle",
                "科学衛星": "https://celestrak.org/NORAD/elements/gp.php?GROUP=science&FORMAT=tle"
            }
            
            selected_url = st.selectbox("衛星グループ", list(url_options.keys()))
            url = url_options[selected_url]
            
            if st.button("TLEデータを取得"):
                with st.spinner("データを取得中..."):
                    tle_data = fetch_tle_data(url)
                    if tle_data:
                        st.session_state.tle_data = tle_data
                        st.session_state.data_loaded = True
                        st.success("データが正常に読み込まれました")
        
        # 手動入力の場合
        elif data_source == "手動入力":
            st.subheader("TLEデータ入力")
            tle_text = st.text_area(
                "TLEデータを入力してください",
                height=200,
                placeholder="衛星名\n1 25544U 98067A   21001.50000000  .00001448  00000-0  33508-4 0  9991\n2 25544  51.6435 114.6994 0003448 202.4638 157.5534 15.48905396276519",
                help="衛星名、TLE行1、TLE行2の順で3行ずつ入力してください"
            )
            
            if st.button("TLEデータを解析"):
                if tle_text.strip():
                    tle_data = parse_tle_text(tle_text)
                    if tle_data:
                        st.session_state.tle_data = tle_data
                        st.session_state.data_loaded = True
                        st.success(f"{len(tle_data)}個の衛星データを解析しました")
                    else:
                        st.error("有効なTLEデータが見つかりませんでした")
                else:
                    st.error("TLEデータを入力してください")
        
        # ローカルファイルの場合
        elif data_source == "TELファイル":
            st.subheader("TELファイル読み込み")
            group_options = {
                "過去30日間": "./tle_data/tle_last30days.txt",
                "気象衛星": "./tle_data/tle_weather.txt",
                "科学衛星": "./tle_data/tle_science.txt"
            }
            selected_group = st.selectbox("衛星グループ", list(group_options.keys()))
            file_path = group_options[selected_group]
            if st.button("TLEデータを読み込む"):
                tle_data = load_tle_from_file(file_path)
                if tle_data:
                    st.session_state.tle_data = tle_data
                    st.session_state.data_loaded = True
                    st.success(f"{len(tle_data)}個の衛星データを読み込みました")
                else:
                    st.error("有効なTLEデータが見つかりませんでした")
        
        # 計算パラメータ
        st.subheader("計算パラメータ")
        end_time = st.number_input(
            "計算点数（120秒間隔）",
            min_value=10,
            max_value=1000,
            value=100,
            step=10,
            help="軌道計算の点数を設定してください(TLEでの時間から何秒後まで計算するか)"
        )
        
        sat_num = st.number_input(
            "表示衛星数",
            min_value=1,
            max_value=100,
            value=20,
            step=1,
            help="表示する衛星の数を設定してください"
        )
        
        # コントロールボタン
        st.subheader("コントロール")
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("計算開始", type="primary"):
                if hasattr(st.session_state, 'data_loaded') and st.session_state.data_loaded:
                    st.session_state.calculate = True
                else:
                    st.error("先にTLEデータを読み込んでください")
        
        with col2:
            if st.button("リセット"):
                if 'tle_data' in st.session_state:
                    del st.session_state.tle_data
                if 'data_loaded' in st.session_state:
                    del st.session_state.data_loaded
                if 'calculate' in st.session_state:
                    del st.session_state.calculate
                st.success("データがリセットされました")
                st.rerun()
    
    # メインエリア
    if hasattr(st.session_state, 'tle_data') and st.session_state.tle_data:
        tabs = st.tabs(["プロット", "衛星リスト"])
        with tabs[1]:
            st.subheader("衛星名一覧")
            parsed_table = []
            for tle in st.session_state.tle_data:
                orbit_elements = parse_tle(tle)
                if orbit_elements:
                    parsed_table.append({
                        '衛星名': orbit_elements['name'][:29],
                        '軌道長半径a': f"{calculate_semi_major_axis(orbit_elements['mean_motion']):.2f} km",
                        '軌道傾斜角i': f"{orbit_elements['inclination']:.4f}°",
                        '昇交点赤経Ω': f"{orbit_elements['raan']:.4f}°",
                        '離心率e': f"{orbit_elements['eccentricity']:.7f}",
                        '近地点引数ω': f"{orbit_elements['arg_perigee']:.4f}°",
                        '平均近点角M': f"{orbit_elements['mean_anomaly']:.4f}°",
                        '平均運動n': f"{orbit_elements['mean_motion']:.8f}°/s",
                        '周期T': f"{orbit_elements['period']/3600:.4f}h"
                    })
            if parsed_table:
                st.dataframe(pd.DataFrame(parsed_table))
            else:
                st.info("パースできるデータがありません。")
        with tabs[0]:
            if hasattr(st.session_state, 'calculate') and st.session_state.calculate:
                with st.spinner("軌道計算中..."):
                    fig = create_orbit_plot(st.session_state.tle_data, end_time, sat_num)
                    st.plotly_chart(fig, use_container_width=True, height=800)
            else:
                st.info("👈 サイドバーからTLEデータを読み込んで計算を開始してください")
    else:
        st.info("👈 サイドバーからTLEデータを読み込んで計算を開始してください")
        
        # 使用例を表示
        st.markdown("""
        ### TLEデータの入力例
        ```
        ISS (ZARYA)
        1 25544U 98067A   21001.50000000  .00001448  00000-0  33508-4 0  9991
        2 25544  51.6435 114.6994 0003448 202.4638 157.5534 15.48905396276519
        ```
        
        TLEデータ
        <https://celestrak.org/NORAD/elements/>
        
        ### 手順
        1. サイドバーでデータソースを選択
        2. TLEデータを取得または入力
        3. 計算パラメータを設定
        4. 「計算開始」ボタンをクリック
        
        ※ニュートン法で軌道を計算
        """)

if __name__ == "__main__":
    main() 