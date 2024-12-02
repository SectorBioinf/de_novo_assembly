import subprocess
import time
import os
import fnmatch
import shutil
import pandas as pd

# Пути для программ
FASTQC = "fastqc * -t 8" 
trimmomatic_jar = 'trimmomatic'
SPADES = "spades.py"
NOVOPLASTY = "NOVOPlasty"
ABYSS_PE = "abyss-pe"
QUAST = "quast.py"
BUSCO = "busco"

# Путь к файлу, где будет сохраняться информация о времени выполнения
LOG_FILE = '/путь/к/time.log'


# Функция для выполнения команды в терминале + запись времени в log-файл
def run_command(command):
    start_time = time.time() #Эта строка записывает текущее время в переменную start_time
    print(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    elapsed_time = time.time() - start_time # Эта строка вычисляет разницу между текущим временем и start_time, чтобы получить время, затраченное на выполнение команды.
    print(f"Command finished in {elapsed_time:.2f} seconds")
    # Записываем время выполнения в файл
    with open(LOG_FILE, 'a') as log_file:
        log_file.write(f"Command: {command}\nElapsed time: {elapsed_time:.2f} seconds\n\n")
        
    return elapsed_time

def quality_control(fastq_files_path, output_dir):
    original_dir = os.getcwd()  # Сохраняем текущую рабочую директорию
    command = f"cd {fastq_files_path} && {FASTQC} -o {output_dir} && cd -"  # Добавляем cd - для возврата в исходную директорию
    run_command(command)
    os.chdir(original_dir)  # Возвращаемся в исходную директорию после выполнения команды
   #print(f"Current directory after cd: {os.getcwd()}")  # для проверки

# Фильтрация

def trim_reads(fastq_files_path, output_dir, trimmomatic_jar):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    # Перебор всех файлов в указанной директории
    for file_name in os.listdir(fastq_files_path):
        if file_name.endswith("_R1_001.fastq.gz"):
            folder_name = file_name.split("_R1")[0]  # Извлечение части имени файла перед _R1

            input_file_R1 = os.path.join(fastq_files_path, file_name)  # Полный путь к файлу R1
            input_file_R2 = os.path.join(fastq_files_path, file_name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs:  # Проверка, была ли пара файлов уже обработана
                output_file_R1_paired = os.path.join(output_dir, f"{folder_name}_paired_R1.fastq.gz")
                output_file_R1_unpaired = os.path.join(output_dir, f"{folder_name}_unpaired_R1.fastq.gz")
                output_file_R2_paired = os.path.join(output_dir, f"{folder_name}_paired_R2.fastq.gz")
                output_file_R2_unpaired = os.path.join(output_dir, f"{folder_name}_unpaired_R2.fastq.gz")
                
                # Команда для запуска Trimmomatic
                command = f"{trimmomatic_jar} PE -phred33 {input_file_R1} {input_file_R2} {output_file_R1_paired} {output_file_R1_unpaired} {output_file_R2_paired} {output_file_R2_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:20 CROP:265"
                
                elapsed_time = run_command(command)  # Запуск команды и получение времени выполнения
                total_elapsed_time += elapsed_time  # Суммирование общего времени
                processed_pairs.add(file_pair)  # Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds")  # Вывод общего времени обработки (1327.12 сек)


# Сборка
# SPAdes
def assemble_spades(input_files, output_dir, kmer_size):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            output_folder = os.path.join(output_dir, folder_name) # Создание пути к выходной папке
            os.makedirs(output_folder, exist_ok=True)  # Создание папки для каждой пары файлов

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = f"{SPADES} -1 {input_file_R1} -2 {input_file_R2} -o {output_folder} -k {kmer_size}"
                elapsed_time = run_command(command)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") #7039.27
    return total_elapsed_time
    
# NOVOPlasty
def assemble_novo(input_folder, output_dir, kmer_size, length, config_template, seed_file_mito, seed_file_chloro):
    total_elapsed_time = 0 

    # Перебор файлов в указанной папке
    for file_name in os.listdir(input_folder):
        if file_name.endswith("_paired_R1.fastq.gz"):
            # Извлечение части имени файла перед _paired
            folder_name = file_name.split("_paired")[0]
            output_folder = os.path.join(output_dir, folder_name) + '/'
            
            # Создание папки для каждой пары файлов
            os.makedirs(output_folder, exist_ok=True)
            
            # Создание путей к файлам
            input_file_R1 = os.path.join(input_folder, file_name)
            input_file_R2 = os.path.join(input_folder, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            
            # Путь к скопированному файлу конфигурации
            config_file_path = os.path.join(output_folder, 'config.txt')
            
            # Копируем шаблонный файл конфигурации в выходную папку
            shutil.copy(config_template, config_file_path)

            # Обработка для митохондрий
            config_type = 'mito_plant'
            genome_range = '400000-500000'
            
            # Обновляем конфигурационный файл для митохондрий
            update_config_file(config_file_path, folder_name, config_type, genome_range, 
                               kmer_size, seed_file_mito, length, input_file_R1, input_file_R2, output_folder)

            # Формируем команду для запуска Novoplasty
            command = f"{NOVOPLASTY} -c {config_file_path}"
            elapsed_time = run_command(command)
            total_elapsed_time += elapsed_time
            
            # Обработка для хлоропластов
            config_type = 'chloro'
            genome_range = '120000-200000'
            
            # Обновляем конфигурационный файл для хлоропластов
            update_config_file(config_file_path, folder_name, config_type, genome_range, 
                               kmer_size, seed_file_chloro, length, input_file_R1, input_file_R2, output_folder)

            # Формируем команду для повторного запуска Novoplasty с новой конфигурацией
            command = f"{NOVOPLASTY} -c {config_file_path}"
            elapsed_time = run_command(command)
            total_elapsed_time += elapsed_time

    print(f"Total elapsed time: {total_elapsed_time} seconds")
    return total_elapsed_time
    
def update_config_file(config_file_path, folder_name, config_type, genome_range, 
                       kmer_size, seed_file, length, input_file_R1, input_file_R2, output_folder):
    # Читаем шаблонный файл конфигурации
    with open(config_file_path, 'r') as file:
        lines = file.readlines()
    
    # Открываем файл для записи
    with open(config_file_path, 'w') as file:
        # Определяем окончание в зависимости от типа конфигурации
        project_suffix = '_mito' if config_type == 'mito_plant' else '_chloro'
        full_project_name = f"{folder_name}{project_suffix}"
        
        for line in lines:
            if line.startswith("Project name"):
                file.write(f"Project name          = {full_project_name}\n")
            elif line.startswith("Type"):
                file.write(f"Type                  = {config_type}\n")
            elif line.startswith("Genome Range"):
                file.write(f"Genome Range          = {genome_range}\n")
            elif line.startswith("K-mer"):
                file.write(f"K-mer                 = {kmer_size}\n")
            elif line.startswith("Seed Input"):
                file.write(f"Seed Input            = {seed_file}\n")
            elif line.startswith("Read Length"):
                file.write(f"Read Length           = {length}\n")
            elif line.startswith("Forward reads"):
                file.write(f"Forward reads         = {input_file_R1}\n")
            elif line.startswith("Reverse reads"):
                file.write(f"Reverse reads         = {input_file_R2}\n")
            elif line.startswith("Output path"):
                file.write(f"Output path           = {output_folder}\n")
            elif line.startswith("Reference sequence") or \
                 line.startswith("Chloroplast sequence") or \
                 line.startswith("Insert size"):
                file.write(f"{line.split('=')[0]} =\n")  # Оставляем только ключевое слово и знак равенства
            else:
                file.write(line)  # Не меняем остальные строки



# ABySS
def assemble_abyss(input_files, output_dir, kmer_size):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            output_folder = os.path.join(output_dir, folder_name) # Создание пути к выходной папке
            os.makedirs(output_folder, exist_ok=True)  # Создание папки для каждой пары файлов

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = f"{ABYSS_PE} k={kmer_size} name={folder_name} in='{input_file_R1} {input_file_R2}' -C {output_folder} j=12 B=30G"
                elapsed_time = run_command(command)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds")
    return total_elapsed_time


#quast

def quast(contigs_dir, output_dir):
    all_contigs_files = []
    
    # Проходим по всем папкам в contigs_dir
    for folder_name in os.listdir(contigs_dir):
        folder_path = os.path.join(contigs_dir, folder_name)
        
        # Проверяем, что это директория
        if os.path.isdir(folder_path):
            # Ищем все файлы, включающие в себя 'contigs.fa'
            for file in os.listdir(folder_path):
                if fnmatch.fnmatch(file, '*contigs.fa*'):
                    all_contigs_files.append(os.path.join(folder_path, file))
    
    # Проверяем, найдены ли файлы
    if all_contigs_files:
        # Создаем основную директорию для вывода, если она не существует
        os.makedirs(output_dir, exist_ok=True)
        
        # Формируем команду с учетом всех файлов
        contigs_files_str = ' '.join(all_contigs_files)
        command = f"{QUAST} {contigs_files_str} -o {output_dir} -t 8"
        
        # Запускаем команду и измеряем время выполнения
        elapsed_time = run_command(command)
        return elapsed_time
    else:
        print("No files matching 'contigs.fa' found in the specified directory.")

# Модификация файла NOVOPlasty для обработки программой QUAST 
def modify_contig_files(contig_files, output_dir):
    modified_contigs = []
    for contig_file in contig_files:
        modified_file = os.path.join(output_dir, os.path.basename(contig_file))
        shutil.copy(contig_file, modified_file)

        # Заменяем '*' на 'N'
        with open(modified_file, 'r') as f:
            content = f.read()
        
        content = content.replace('*', 'N')

        with open(modified_file, 'w') as f:
            f.write(content)
        
        # Добавляем модифицированный файл в список для команды
        modified_contigs.append(modified_file)
    
    return modified_contigs

def quast_novo(contigs_dir, output_dir):
    contigs_to_modify_mito = []
    contigs_to_modify_chloro = []

    # Проходим по всем папкам в contigs_dir
    for folder_name in os.listdir(contigs_dir):
        folder_path = os.path.join(contigs_dir, folder_name)

        # Проверяем, что это директория
        if os.path.isdir(folder_path):
            # Первый цикл для митохондриальных файлов
            found_circularized_mito = False

            for file in os.listdir(folder_path):
                file_path = os.path.join(folder_path, file)
    
                if fnmatch.fnmatch(file, 'Circularized*_mito.fasta'):
                    contigs_to_modify_mito.append(file_path)
                    found_circularized_mito = True

            # Проверяем наличие Contigs, только если файлы Circularized не найдены
            if not found_circularized_mito:
                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
                    if fnmatch.fnmatch(file, 'Contigs*_mito.fasta'):
                        contigs_to_modify_mito.append(file_path)

            # Второй цикл для хлоропластных файлов
            found_circularized_chloro = False

            for file in os.listdir(folder_path):
                file_path = os.path.join(folder_path, file)
    
                if fnmatch.fnmatch(file, 'Circularized*_chloro.fasta'):
                    contigs_to_modify_chloro.append(file_path)
                    found_circularized_chloro = True

            # Проверяем наличие Contigs, только если файлы Circularized не найдены
            if not found_circularized_chloro:
                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
                    if fnmatch.fnmatch(file, 'Contigs*_chloro.fasta'):
                        contigs_to_modify_chloro.append(file_path)
                    

    # Проверка на наличие файлов
    if not contigs_to_modify_mito and not contigs_to_modify_chloro:
        print("No files matching 'Contigs*.fasta' or 'Circularized*.fasta' found in the specified directory.")
        return

    # Создаем основную директорию для вывода, если она не существует
    mito_output_dir = os.path.join(output_dir, 'mito')
    chloro_output_dir = os.path.join(output_dir, 'chloro')
    os.makedirs(mito_output_dir, exist_ok=True)
    os.makedirs(chloro_output_dir, exist_ok=True)

    # Модифицируем файлы и получаем списки модифицированных файлов
    modified_contigs_mito = modify_contig_files(contigs_to_modify_mito, mito_output_dir)
    modified_contigs_chloro = modify_contig_files(contigs_to_modify_chloro, chloro_output_dir)

    # Формируем команды для _mito и _chloro
    command_mito = f"{QUAST} {' '.join(modified_contigs_mito)} -o {mito_output_dir} -t 12"
    command_chloro = f"{QUAST} {' '.join(modified_contigs_chloro)} -o {chloro_output_dir} -t 12"

    # Запускаем команды и измеряем время выполнения
    total_time = run_command(command_mito) + run_command(command_chloro)
    print (f"Total time for processing: {total_time:.2f} seconds")
    return total_time

#busco

def busco(scaffolds_dir, output_dir, lineage_dataset, code_dir):
    all_scaffolds_files = []
    temp_scaffolds_dir = os.path.join(output_dir, "temp_scaffolds")
    result_dir = os.path.join(output_dir, "result")

    # Проходим по всем папкам в scaffolds_dir
    for folder_name in os.listdir(scaffolds_dir):
        folder_path = os.path.join(scaffolds_dir, folder_name)

        # Проверяем, что это директория
        if os.path.isdir(folder_path):
            # Ищем файлы с именем 'scaffolds.fa' и 'Contigs*.fasta'
            for file in os.listdir(folder_path):
                if fnmatch.fnmatch(file, '*scaffolds.fa*') or fnmatch.fnmatch(file, 'Contigs*.fasta'):
                    # Добавляем найденные файлы в список
                    all_scaffolds_files.append(os.path.join(folder_path, file))
                    
    print (all_scaffolds_files)
    # Проверяем, найдены ли файлы
    if all_scaffolds_files:
        # Создаем временную директорию для копирования файлов, если она не существует
        os.makedirs(temp_scaffolds_dir, exist_ok=True)

        # Копируем все найденные файлы во временную директорию с изменением имени
        for file in all_scaffolds_files:
            # Извлекаем имя родительской папки для использования в новом имени файла
            folder_name = os.path.basename(os.path.dirname(file))
            # Получаем только расширение файла
            _, extension = os.path.splitext(os.path.basename(file))
            # Формируем новое имя файла с расширением
            new_file_name = f"{folder_name}{extension}"
            # Копируем файл во временную директорию с новым именем
            shutil.copy(file, os.path.join(temp_scaffolds_dir, new_file_name))     

        # Формируем команду с учетом временной директории
        command = f"{BUSCO} -i {temp_scaffolds_dir} -c 12 --out_path {output_dir} --mode genome --offline -l {lineage_dataset} --metaeuk"
        
        # Запускаем команду и измеряем время выполнения
        elapsed_time_busco = run_command(command)
        
        # Удаляем временную директорию после завершения анализа
        shutil.rmtree(temp_scaffolds_dir)
        
        # Ищем файлы результатов
        all_summary_files = []
        
        # Проходим по всем папкам в output_dir
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if fnmatch.fnmatch(file, 'short_summary.*.txt'):
                    all_summary_files.append(os.path.join(root, file))
        
        # Проверяем, найдены ли файлы
        if all_summary_files:
            # Создаем директорию для копирования файлов
            os.makedirs(result_dir, exist_ok=True)
            
            # Копируем все найденные файлы в директорию
            for file in all_summary_files:
                shutil.copy(file, result_dir)
            
            # Выводим содержимое директории
            print("Содержимое директории result:")
            print(os.listdir(result_dir))
            
            # Формируем команду для запуска
            command_plot = f"python3 {code_dir} -wd {result_dir}"
            
            # Запускаем команду и измеряем время выполнения
            elapsed_time_plot = run_command(command_plot)
            
            # Объединяем время
            total_elapsed_time = elapsed_time_busco + elapsed_time_plot
            print(f"Total elapsed time for BUSCO and plot generation: {total_elapsed_time:.2f} seconds.")
            return total_elapsed_time

#Модификация файлов NOVOPlasty для обработки программой BUSCO
def process_file(file_path, temp_scaffolds_dir):
    # Открываем текущий файл для чтения
    with open(file_path, 'r') as f:
        lines = f.readlines()  # Читаем все строки файла

    # Словарь для отслеживания количества заголовков
    seen_counts = {}
    new_lines = []

    # Проходим по каждой строке в файле
    for line in lines:
        # Заменяем '*' на 'N' в строке
        line = line.replace('*', 'N')
        if line.startswith('>'):  # Проверяем, начинается ли строка с символа '>'
            if line in seen_counts:
                seen_counts[line] += 1  # Увеличиваем счетчик для дубликата
                new_line = f"{line.strip()}({seen_counts[line]})\n"
                new_lines.append(new_line)  # Добавляем обновленную строку в список
            else:
                seen_counts[line] = 1  # Если уникальная, добавляем ее в словарь
                new_lines.append(line)  # Добавляем строку без изменений
        else:
            new_lines.append(line)  # Если не заголовок, добавляем без изменений

    # Сохраняем обновленный файл во временной директории
    output_file_path = os.path.join(temp_scaffolds_dir, os.path.basename(file_path))
    with open(output_file_path, 'w') as f_out:
        f_out.writelines(new_lines)  # Записываем все обновленные строки в новый файл 
           
def busco_novo(scaffolds_dir, output_dir, lineage_dataset, code_dir):
    temp_scaffolds_dir = os.path.join(output_dir, "temp_scaffolds")
    result_dir = os.path.join(output_dir, "result")
    
    # Создаем временную директорию, если она не существует
    os.makedirs(temp_scaffolds_dir, exist_ok=True)

    # Проходим по всем папкам в scaffolds_dir
    for folder_name in os.listdir(scaffolds_dir):
        folder_path = os.path.join(scaffolds_dir, folder_name)

        # Проверяем, что это директория
        if os.path.isdir(folder_path):
            contigs_to_modify_mito = []
            contigs_to_modify_chloro = []

            # Проверяем, что это директория
            if os.path.isdir(folder_path):
                # Первый цикл для митохондриальных файлов
                found_circularized_mito = False

                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
    
                    if fnmatch.fnmatch(file, 'Circularized*_mito.fasta'):
                        contigs_to_modify_mito.append(file_path)
                        found_circularized_mito = True

                # Проверяем наличие Contigs, только если файлы Circularized не найдены
                if not found_circularized_mito:
                    for file in os.listdir(folder_path):
                        file_path = os.path.join(folder_path, file)
                        if fnmatch.fnmatch(file, 'Contigs*_mito.fasta'):
                            contigs_to_modify_mito.append(file_path)

                # Второй цикл для хлоропластных файлов
                found_circularized_chloro = False

                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
    
                    if fnmatch.fnmatch(file, 'Circularized*_chloro.fasta'):
                        contigs_to_modify_chloro.append(file_path)
                        found_circularized_chloro = True

                # Проверяем наличие Contigs, только если файлы Circularized не найдены
                if not found_circularized_chloro:
                    for file in os.listdir(folder_path):
                        file_path = os.path.join(folder_path, file)
                        if fnmatch.fnmatch(file, 'Contigs*_chloro.fasta'):
                            contigs_to_modify_chloro.append(file_path)
   
            # Обработка файлов для митохондрий
            for file_path in contigs_to_modify_mito:
                process_file(file_path, temp_scaffolds_dir)

            # Обработка файлов для хлоропластов
            for file_path in contigs_to_modify_chloro:
                process_file(file_path, temp_scaffolds_dir)

    # Формируем команду с учетом временной директории
    command = f"{BUSCO} -i {temp_scaffolds_dir} -c 12 --out_path {output_dir} --mode genome --offline -l {lineage_dataset} --metaeuk"
    elapsed_time_busco = run_command(command)
    shutil.rmtree(temp_scaffolds_dir)
    mito_files = []
    chloro_files = []
    
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if fnmatch.fnmatch(file, 'short_summary.*mito*.txt'):
                mito_files.append(os.path.join(root, file))
            elif fnmatch.fnmatch(file, 'short_summary.*chloro*.txt'):
                chloro_files.append(os.path.join(root, file))
    
     # Создание директории для копирования
    if mito_files or chloro_files:
        os.makedirs(os.path.join(result_dir, 'mito'), exist_ok=True)
        os.makedirs(os.path.join(result_dir, 'chloro'), exist_ok=True)

        # Копирование файлов в директорию mito
        for file in mito_files:
            shutil.copy(file, os.path.join(result_dir, 'mito'))

        # Копирование файлов в директорию chloro
        for file in chloro_files:
            shutil.copy(file, os.path.join(result_dir, 'chloro'))
        
        # Выводим содержимое директории
        print("Содержимое директории result:")
        print(os.listdir(result_dir))
        
        command_mito = f"python3 {code_dir} -wd {os.path.join(result_dir, 'mito')}"
        command_chloro = f"python3 {code_dir} -wd {os.path.join(result_dir, 'chloro')}"
        elapsed_time_plot = run_command(command_mito)+run_command(command_chloro)
        total_elapsed_time = elapsed_time_busco + elapsed_time_plot
        print(f"Общее время выполнения BUSCO и генерации графиков: {total_elapsed_time:.2f} секунд.")
        return total_elapsed_time

    
# Функция для создания отчета в excel
def create_excel_report(output_file, program_info):
    # Преобразуем данные в DataFrame
    df = pd.DataFrame(list(program_info.items()), columns=['Название программы', 'Время выполнения (с)'])

    # Проверяем, существует ли файл
    if os.path.exists(output_file):
        # Если файл существует, то добавляем данные в конец
        with pd.ExcelWriter(output_file, engine='openpyxl', mode='a', if_sheet_exists='overlay') as writer:
            # Найдем последнюю строку, чтобы добавить новые данные после старых
            startrow = writer.sheets['Sheet1'].max_row
            df.to_excel(writer, index=False, header=False, startrow=startrow)
        print(f"Новые данные добавлены в {output_file}")
    else:
        # Если файл не существует, создаем новый и записываем данные
        df.to_excel(output_file, index=False)
        print(f"Отчет сохранен в {output_file}")
    
# Main
def main():
    # Пути к входным и выходным данным
    fastq_files_path = '/путь/к/Fastq'
    output_dir = '/путь/к/Fastq_output'
    #trimmed_files = '/путь/к/Trim'
    
    trimmed_files_good = '/путь/к/good_reads'
    trimmed_files_bad = '/путь/к/bad_reads'
    
    novo_output_good = '/путь/к/novo_good'
    novo_output_bad = '/путь/к/novo_bad'
    config_template = '/путь/к/config.txt'
    
        # Путь к seed файлам для NOVOPlasty
    seed_file_mito = '/путь/к/seed.fasta'
    seed_file_chloro = '/путь/к/Seed_RUBP.fasta'
    
    spades_output_good = '/путь/к/spades_good'
    spades_output_bad = '/путь/к/spades_bad'
    
    abyss_output_good = '/путь/к/abyss_good'
    abyss_output_bad = '/путь/к/abyss_bad'
    
    quast_output_spades_good = '/путь/к/spades_good'
    quast_output_spades_bad = '/путь/к/spades_bad'
    
    quast_output_abyss_good = '/путь/к/abyss_good'
    quast_output_abyss_bad = '/путь/к/abyss_bad'
    
    quast_output_novo_good = '/путь/к/novo_good'
    quast_output_novo_bad = '/путь/к/novo_bad'
    
      #Пути к библиотеке для BUSCO и к коду для построения диаграммы 
    lineage_dataset = '/путь/к/lineage_dataset'
    code_dir = '/путь/к/generate_plot.py'
    
    busco_output_spades_good = '/путь/к/spades_good'
    busco_output_spades_bad = '/путь/к/spades_bad'
    
    busco_output_abyss_good = '/путь/к/abyss_good'
    busco_output_abyss_bad = '/путь/к/abyss_bad'
    
    busco_output_novo_good = '/путь/к/novo_good2'
    busco_output_novo_bad = '/путь/к/novo_bad'
    
    
    # Запуск функций
    # QC
    #quality_control(fastq_files_path, output_dir)
    
    # Trimming
    #trim_reads(fastq_files_path, trimmed_files, trimmomatic_jar)
    
    # Assembling 
    #spades_good = assemble_spades(trimmed_files_good, spades_output_good, 55)
    #spades_bad = assemble_spades(trimmed_files_bad, spades_output_bad, 55)
    
    #novo_good = assemble_novo(trimmed_files_good, novo_output_good, 55, 265, config_template, seed_file_mito, seed_file_chloro)
    #novo_bad = assemble_novo(trimmed_files_bad, novo_output_bad, 55, 265, config_template, seed_file_mito, seed_file_chloro)
    
    #abyss_good = assemble_abyss(trimmed_files_good, abyss_output_good, 64)
    #abyss_bad = assemble_abyss(trimmed_files_bad, abyss_output_bad, 64)

    #quast_spades_good = quast(spades_output_good, quast_output_spades_good)
    #quast_spades_bad = quast(spades_output_bad, quast_output_spades_bad)
    
    #quast_abyss_good = quast(abyss_output_good, quast_output_abyss_good)
    #quast_abyss_bad = (abyss_output_bad, quast_output_abyss_bad)
    
    #quast_novo_good = quast_novo(novo_output_good, quast_output_novo_good)
    #quast_novo_bad = quast_novo(novo_output_bad, quast_output_novo_bad)
    
    #busco_spades_good = busco(spades_output_good, busco_output_spades_good, lineage_dataset,code_dir)
    #busco_spades_bad = busco(spades_output_bad, busco_output_spades_bad, lineage_dataset,code_dir)
    
    #busco_abyss_good = busco(abyss_output_good, busco_output_abyss_good, lineage_dataset,code_dir)
    #busco_abyss_bad = busco(abyss_output_bad, busco_output_abyss_bad, lineage_dataset,code_dir)
    
    #busco_novo_good = busco_novo(novo_output_good, busco_output_novo_good, lineage_dataset,code_dir)
    #busco_novo_bad = busco_novo(novo_output_bad, busco_output_novo_bad, lineage_dataset,code_dir)
    
    #создание таблицы excel
    #program_info = {}
    #program_info['SPAdes_good'] = spades_good
    #program_info['SPAdes_bad'] = spades_bad
    #program_info['NOVOPlasty_good'] = novo_good
    #program_info['NOVOPlasty_bad'] = novo_bad
    #program_info['ABySS_good'] = abyss_good
    #program_info['ABySS_bad'] = abyss_good
    #program_info['QUAST_spades_good'] = quast_spades_good
    #program_info['QUAST_spades_bad'] = quast_spades_bad
    #program_info['QUAST_abyss_good'] = quast_abyss_good
    #program_info['QUAST_abyss_bad'] = quast_abyss_bad
    #program_info['QUAST_novo_good'] = quast_novo_good
    #program_info['QUAST_novo_bad'] = quast_novo_bad
    #program_info['BUSCO_spades_good'] = busco_spades_good
    #program_info['BUSCO_spades_bad'] = busco_spades_bad
    #program_info['BUSCO_abyss_good'] = busco_abyss_good
    #program_info['BUSCO_abyss_bad'] = busco_abyss_bad
    #program_info['BUSCO_novo_good'] = busco_novo_good
    #program_info['BUSCO_novo_bad'] = busco_novo_bad
    
    # Создаем отчет
    #create_excel_report('quast.xlsx', program_info)

if __name__ == "__main__":
    main()



