# Используем официальный Python образ
FROM python:3.9

# Устанавливаем рабочую директорию
WORKDIR /app

# Копируем файлы
COPY requirements.txt .
COPY dash_app.py .
COPY streamlit_app.py .
COPY supporting_files .

# Устанавливаем зависимости
RUN pip install --no-cache-dir -r requirements.txt

# Открываем порты
EXPOSE 8501 8700

# Запускаем оба сервера
CMD python dash_app.py & streamlit run streamlit_app.py --server.port=8501 --server.address=0.0.0.0
