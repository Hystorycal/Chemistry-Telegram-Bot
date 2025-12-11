from dotenv import load_dotenv # Импорт библиотеки
load_dotenv() # Загрузка переменных из .env
import logging
import os
import requests
import matplotlib

# Устанавливаем бекенд 'Agg' (для работы без окна)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import io
import numpy as np

# Химия
from chempy import Substance, balance_stoichiometry
import mendeleev
import pubchempy as pcp

from telegram import Update, ReplyKeyboardMarkup, ReplyKeyboardRemove
from telegram.ext import (
    Application,

    CommandHandler,
    MessageHandler,
    filters,
    ContextTypes,
    ConversationHandler,
)

# Логирование
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# Состояния
CHOOSING, MOLAR_MASS, EQUALIZE, VALENCY, INDICATORS, RESIDUES, CONFIG, HYBRIDIZATION = range(8)

# Кнопки
BTN_MOLAR = "Calculate molar mass and give image of structure"
BTN_EQUALIZE = "Equalize the reaction"
BTN_VALENCY = "Write all possible oxidation states"
BTN_INDICATORS = "Indicator Colors (Lac/Phen/Meth)"
BTN_RESIDUES = "Residues (Color/State)"
BTN_CONFIG = "Electronic configuration"
BTN_HYBRID = "Hybridization"

MAIN_KEYBOARD = [
    [BTN_MOLAR, BTN_EQUALIZE],
    [BTN_VALENCY, BTN_CONFIG],
    [BTN_INDICATORS, BTN_RESIDUES],
    [BTN_HYBRID]
]


# --- 1. ГЕНЕРАТОР ОРБИТАЛЕЙ (ВСТРОЕН В КОД) ---
def create_lobe(ax, angle_deg, color='skyblue'):
    """Рисует один лепесток орбитали под нужным углом"""
    angle_rad = np.radians(angle_deg)

    # Смещаем центр эллипса, чтобы он исходил из центра атома
    distance = 0.6
    x = distance * np.cos(angle_rad)
    y = distance * np.sin(angle_rad)

    # Рисуем эллипс (лепесток)
    ell = Ellipse((x, y), width=1.2, height=0.5, angle=angle_deg,
                  facecolor=color, edgecolor='black', alpha=0.8, zorder=2)
    ax.add_patch(ell)

    # Маленький "хвостик" (противофаза) с другой стороны
    x_tail = -0.2 * np.cos(angle_rad)
    y_tail = -0.2 * np.sin(angle_rad)
    ell_tail = Ellipse((x_tail, y_tail), width=0.4, height=0.3, angle=angle_deg,
                       facecolor='white', edgecolor='gray', alpha=0.5, zorder=1)
    ax.add_patch(ell_tail)


def generate_orbital_images():
    """Создает картинки гибридизации (sp, sp2...) локально"""
    print("🎨 Generating orbital diagrams...")

    # Конфигурация углов для разных типов
    configs = {
        "sp": [0, 180],
        "sp2": [90, 210, 330],
        "sp3": [90, 200, 340, 270],  # 2D проекция тетраэдра
        "sp3d": [90, 270, 0, 120, 240],  # Тригональная бипирамида (упрощенно)
        "sp3d2": [0, 60, 120, 180, 240, 300]  # Октаэдр (проекция)
    }

    colors = ['#FFCC00', '#00CCFF', '#FF6666', '#99CC00', '#CC99FF', '#FF99CC']

    for name, angles in configs.items():
        filename = f"{name}.png"
        if os.path.exists(filename):
            continue

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_aspect('equal')
        ax.axis('off')

        # Рисуем лепестки
        for i, angle in enumerate(angles):
            color = colors[i % len(colors)]
            create_lobe(ax, angle, color)

        # Рисуем ядро
        nucleus = plt.Circle((0, 0), 0.15, color='black', zorder=10)
        ax.add_patch(nucleus)

        # Подпись
        ax.text(0, -1.8, f"{name} hybridization", ha='center', fontsize=16, fontweight='bold')

        # Настройка границ
        ax.set_xlim(-2, 2)
        ax.set_ylim(-2, 2)

        plt.tight_layout()
        fig.savefig(filename, dpi=100, bbox_inches='tight')
        plt.close(fig)
        print(f"✅ Generated {filename}")


# --- ОСТАЛЬНЫЕ ФУНКЦИИ ---

def get_pubchem_2d_image(query):
    try:
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query.strip()}/cids/JSON"
        response = requests.get(search_url, timeout=5)
        if response.status_code != 200:
            search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/{query.strip()}/cids/JSON"
            response = requests.get(search_url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            cid = data['IdentifierList']['CID'][0]
            return f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?record_type=2d&image_size=large"
    except:
        pass
    return None


def draw_atom_structure(symbol, name, atomic_number):
    try:
        electrons_remaining = atomic_number
        shells = []
        n = 1
        while electrons_remaining > 0:
            if n == 1:
                cap = 2
            elif n == 2:
                cap = 8
            elif n == 3 and atomic_number <= 20:
                cap = 8
            else:
                cap = 2 * n ** 2
            count = min(cap, electrons_remaining)
            shells.append(count)
            electrons_remaining -= count
            n += 1
    except:
        shells = [atomic_number]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.axis('off')
    nucleus = plt.Circle((0, 0), 0.5, color='orange', label='Nucleus')
    ax.add_patch(nucleus)
    ax.text(0, 0, f"{symbol}\nZ={atomic_number}", ha='center', va='center', fontweight='bold')
    max_radius = len(shells) + 1
    ax.set_xlim(-max_radius - 1, max_radius + 1)
    ax.set_ylim(-max_radius - 1, max_radius + 1)
    for i, electron_count in enumerate(shells):
        radius = i + 1.5
        circle = plt.Circle((0, 0), radius, fill=False, color='gray', linestyle='--')
        ax.add_patch(circle)
        angles = np.linspace(0, 2 * np.pi, electron_count, endpoint=False)
        for angle in angles:
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            electron = plt.Circle((x, y), 0.15, color='blue')
            ax.add_patch(electron)

    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    plt.close(fig)
    return buf


# --- HANDLERS ---

async def show_menu(update: Update, text: str):
    await update.message.reply_text(
        text,
        reply_markup=ReplyKeyboardMarkup(MAIN_KEYBOARD, one_time_keyboard=True, input_field_placeholder="Select option")
    )
    return CHOOSING


async def start(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    return await show_menu(update, "Hello! I am your Chemistry Assistant.\n\nWhat do you want me to do?")


async def route_selection(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    user_choice = update.message.text
    context.user_data['choice'] = user_choice

    if user_choice == BTN_MOLAR:
        await update.message.reply_text("Please enter chemical formula or name (e.g., H2SO4, Glucose):")
        return MOLAR_MASS
    elif user_choice == BTN_EQUALIZE:
        await update.message.reply_text("Enter equation (e.g., H2 + O2 = H2O):")
        return EQUALIZE
    elif user_choice == BTN_VALENCY:
        await update.message.reply_text("Enter Element Symbol (e.g., Fe, Cl):")
        return VALENCY
    elif user_choice == BTN_INDICATORS:
        await update.message.reply_text("Enter medium (e.g., Acid, Base):")
        return INDICATORS
    elif user_choice == BTN_RESIDUES:
        await update.message.reply_text("Enter precipitate formula (e.g., AgCl):")
        return RESIDUES
    elif user_choice == BTN_CONFIG:
        await update.message.reply_text("Enter Element Symbol (e.g., Ca, Au):")
        return CONFIG
    elif user_choice == BTN_HYBRID:
        await update.message.reply_text("Enter molecule formula (e.g., CH4, BeCl2):")
        return HYBRIDIZATION
    else:
        return await show_menu(update, "Please select a button.")


async def handle_hybridization(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    formula = update.message.text.strip()

    if formula in [BTN_MOLAR, BTN_EQUALIZE, BTN_VALENCY, BTN_INDICATORS, BTN_RESIDUES, BTN_CONFIG, BTN_HYBRID]:
        await update.message.reply_text("⚠️ Please type the formula (e.g., CH4), don't press the menu button.")
        return HYBRIDIZATION

    data_map = {
        "CO2": ("sp", "Linear"), "BeCl2": ("sp", "Linear"), "C2H2": ("sp", "Linear"),
        "BF3": ("sp2", "Trigonal Planar"), "SO3": ("sp2", "Trigonal Planar"), "C2H4": ("sp2", "Trigonal Planar"),
        "CH4": ("sp3", "Tetrahedral"), "NH3": ("sp3", "Pyramidal"), "H2O": ("sp3", "Bent"),
        "PCl5": ("sp3d", "Trigonal Bipyramidal"), "SF6": ("sp3d2", "Octahedral")
    }

    info = data_map.get(formula)

    if info:
        hyb_type = info[0]
        desc = info[1]

        await update.message.reply_text(
            f"ℹ️ **Molecule**: {formula}\n"
            f"⚛️ **Hybridization**: {hyb_type}\n"
            f"📐 **Geometry**: {desc}",
            parse_mode="Markdown"
        )

        # Отправляем локальный файл
        filename = f"{hyb_type}.png"

        if os.path.exists(filename):
            try:
                with open(filename, 'rb') as photo_file:
                    await update.message.reply_photo(photo=photo_file)
            except Exception as e:
                await update.message.reply_text(f"⚠️ Error sending file: {e}")
        else:
            await update.message.reply_text("⚠️ Diagram file generation failed.")

    else:
        await update.message.reply_text(f"⚠️ Unknown molecule: '{formula}'. Try CH4, NH3, BeCl2...")

    return await show_menu(update, "What else do you want to do?")


# --- ОСТАЛЬНЫЕ (Molar, Equalize...) ---

async def handle_molar_mass(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    user_input = update.message.text.strip()
    mass = None;
    formula = user_input;
    source = "Calculation"
    try:
        compounds = pcp.get_compounds(user_input, 'name')
        if not compounds: compounds = pcp.get_compounds(user_input, 'formula')
        if compounds:
            c = compounds[0];
            mass = float(c.molecular_weight);
            formula = c.molecular_formula;
            source = "PubChem"
    except:
        pass
    if mass is None:
        try:
            sub = Substance.from_formula(user_input); mass = sub.mass; source = "ChemPy"
        except:
            pass

    if mass is not None:
        msg = f"🧪 **Input**: {user_input}\n📝 **Formula**: {formula}\n⚖️ **Molar Mass**: {mass:.3f} g/mol\n🔍 **Source**: {source}"
        await update.message.reply_text(msg, parse_mode="Markdown")
        img_url = get_pubchem_2d_image(user_input)
        if img_url: await update.message.reply_photo(photo=img_url)
    else:
        await update.message.reply_text(f"❌ Not found: {user_input}")
    return await show_menu(update, "What else?")


async def handle_equalize(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    text = update.message.text
    try:
        sides = text.split("->") if "->" in text else text.split("=")
        reac = {x.strip() for x in sides[0].split("+") if x.strip()};
        prod = {x.strip() for x in sides[1].split("+") if x.strip()}
        rb, pb = balance_stoichiometry(reac, prod)
        lhs = " + ".join([f"{v if v > 1 else ''}{k}" for k, v in rb.items()]);
        rhs = " + ".join([f"{v if v > 1 else ''}{k}" for k, v in pb.items()])
        await update.message.reply_text(f"✅ Balanced:\n{lhs} → {rhs}")
    except Exception as e:
        await update.message.reply_text(f"❌ Error: {e}")
    return await show_menu(update, "What else?")


async def handle_valency(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    try:
        elem = mendeleev.element(update.message.text.strip())
        ox = elem.oxistates
        await update.message.reply_text(
            f"ℹ️ Oxidation states ({elem.name}):\n{', '.join(map(str, ox)) if ox else 'Unknown'}")
    except:
        await update.message.reply_text("❌ Element not found.")
    return await show_menu(update, "What else?")


async def handle_indicators(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    inp = update.message.text.lower();
    state = "neutral"
    if inp in ["acid", "acidic"] or (inp.startswith('h') and inp != 'h2o'):
        state = "acid"
    elif inp in ["base", "basic"] or "oh" in inp:
        state = "base"
    if state == "acid":
        msg = "🧪 ACID:\n🔴 Litmus: Red\n⚪ Phenolphthalein: Colorless\n🔴 Methyl Orange: Red"
    elif state == "base":
        msg = "🧪 BASE:\n🔵 Litmus: Blue\n🟣 Phenolphthalein: Pink\n🟡 Methyl Orange: Yellow"
    else:
        msg = "🧪 NEUTRAL:\n🟣 Litmus: Violet\n⚪ Phenolphthalein: Colorless\n🟠 Methyl Orange: Orange"
    await update.message.reply_text(msg);
    return await show_menu(update, "What else?")


async def handle_residues(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    db = {"AgCl": "White ppt", "AgBr": "Pale Yellow ppt", "AgI": "Yellow ppt", "BaSO4": "White ppt",
          "PbI2": "Golden Yellow ppt", "Cu(OH)2": "Blue ppt", "Fe(OH)2": "Dirty Green ppt",
          "Fe(OH)3": "Reddish Brown ppt", "CuO": "Black solid"}
    await update.message.reply_text(
        f"ℹ️ {update.message.text}: {db.get(update.message.text.strip(), 'Data not found')}");
    return await show_menu(update, "What else?")


async def handle_config(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    try:
        elem = mendeleev.element(update.message.text.strip())
        await update.message.reply_text(f"⚛️ Config ({elem.name}):\n{elem.econf}")
        await update.message.reply_photo(photo=draw_atom_structure(elem.symbol, elem.name, elem.atomic_number))
    except:
        await update.message.reply_text("❌ Invalid symbol.")
    return await show_menu(update, "What else?")


async def cancel(update: Update, context: ContextTypes.DEFAULT_TYPE) -> int:
    await update.message.reply_text("Cancelled.", reply_markup=ReplyKeyboardRemove());
    return ConversationHandler.END


# --- ЗАПУСК ---

def main() -> None:
    TOKEN = os.getenv("BOT_TOKEN")
    # 1. Сначала сами рисуем картинки (без скачивания!)
    generate_orbital_images()

    # 2. Запускаем бота
    application = Application.builder().token(TOKEN).build()

    conv_handler = ConversationHandler(
        entry_points=[CommandHandler("start", start)],
        states={
            CHOOSING: [MessageHandler(filters.TEXT & ~filters.COMMAND, route_selection)],
            MOLAR_MASS: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_molar_mass)],
            EQUALIZE: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_equalize)],
            VALENCY: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_valency)],
            INDICATORS: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_indicators)],
            RESIDUES: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_residues)],
            CONFIG: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_config)],
            HYBRIDIZATION: [MessageHandler(filters.TEXT & ~filters.COMMAND, handle_hybridization)],
        },
        fallbacks=[CommandHandler("cancel", cancel)],
    )
    application.add_handler(conv_handler)
    print("Bot is running...")
    application.run_polling()


if __name__ == "__main__":
    main()