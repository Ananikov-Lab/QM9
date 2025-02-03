FROM python:3.9

WORKDIR /app

# Копируем зависимости и устанавливаем их
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Копируем весь код в контейнер
COPY . .

# Открываем порты
EXPOSE 5010 5011

# Запуск Streamlit
CMD ["streamlit", "run", "streamlit_app.py", "--server.port=5010", "--server.address=0.0.0.0", "--server.enableCORS=false", "--server.enableXsrfProtection=false"]

