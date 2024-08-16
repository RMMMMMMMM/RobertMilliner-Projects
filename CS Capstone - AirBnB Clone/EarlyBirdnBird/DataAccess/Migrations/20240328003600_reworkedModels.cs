using System;
using Microsoft.EntityFrameworkCore.Migrations;

#nullable disable

namespace EarlyBirdnBird.Data.Migrations
{
    /// <inheritdoc />
    public partial class reworkedModels : Migration
    {
        /// <inheritdoc />
        protected override void Up(MigrationBuilder migrationBuilder)
        {
            migrationBuilder.DropForeignKey(
                name: "FK_Images_Rentals_RentalID",
                table: "Images");

            migrationBuilder.DropForeignKey(
                name: "FK_Rentals_Amenities_AmenityId",
                table: "Rentals");

            migrationBuilder.DropForeignKey(
                name: "FK_Rentals_Beds_BedId",
                table: "Rentals");

            migrationBuilder.DropForeignKey(
                name: "FK_Rentals_Blocked_BlockId",
                table: "Rentals");

            migrationBuilder.DropForeignKey(
                name: "FK_Rentals_Images_ImageId",
                table: "Rentals");

            migrationBuilder.DropForeignKey(
                name: "FK_Rentals_Prices_PriceId",
                table: "Rentals");

            migrationBuilder.DropForeignKey(
                name: "FK_Reservation_Rentals_RentalId",
                table: "Reservation");

            migrationBuilder.DropIndex(
                name: "IX_Reservation_RentalId",
                table: "Reservation");

            migrationBuilder.DropIndex(
                name: "IX_Rentals_AmenityId",
                table: "Rentals");

            migrationBuilder.DropIndex(
                name: "IX_Rentals_BedId",
                table: "Rentals");

            migrationBuilder.DropIndex(
                name: "IX_Rentals_BlockId",
                table: "Rentals");

            migrationBuilder.DropIndex(
                name: "IX_Rentals_ImageId",
                table: "Rentals");

            migrationBuilder.DropIndex(
                name: "IX_Rentals_PriceId",
                table: "Rentals");

            migrationBuilder.DropIndex(
                name: "IX_Images_RentalID",
                table: "Images");

            migrationBuilder.DropColumn(
                name: "AmenityId",
                table: "Rentals");

            migrationBuilder.DropColumn(
                name: "BedId",
                table: "Rentals");

            migrationBuilder.DropColumn(
                name: "BlockId",
                table: "Rentals");

            migrationBuilder.DropColumn(
                name: "ImageId",
                table: "Rentals");

            migrationBuilder.DropColumn(
                name: "Primary",
                table: "Images");

            migrationBuilder.DropColumn(
                name: "Quantity",
                table: "Beds");

            migrationBuilder.RenameColumn(
                name: "RentalId",
                table: "Reservation",
                newName: "Rental");

            migrationBuilder.RenameColumn(
                name: "PriceId",
                table: "Rentals",
                newName: "Beds");

            migrationBuilder.RenameColumn(
                name: "RentalID",
                table: "Images",
                newName: "RentalId");

            migrationBuilder.AlterColumn<double>(
                name: "SalesTax",
                table: "Reservation",
                type: "float",
                nullable: true,
                oldClrType: typeof(float),
                oldType: "real",
                oldNullable: true);

            migrationBuilder.AlterColumn<double>(
                name: "IndividualTotal",
                table: "Reservation",
                type: "float",
                nullable: true,
                oldClrType: typeof(float),
                oldType: "real",
                oldNullable: true);

            migrationBuilder.AlterColumn<double>(
                name: "GrandTotal",
                table: "Reservation",
                type: "float",
                nullable: false,
                oldClrType: typeof(float),
                oldType: "real");

            migrationBuilder.AlterColumn<string>(
                name: "RentalType",
                table: "Rentals",
                type: "nvarchar(max)",
                nullable: false,
                defaultValue: "",
                oldClrType: typeof(string),
                oldType: "nvarchar(max)",
                oldNullable: true);

            migrationBuilder.AddColumn<bool>(
                name: "Published",
                table: "Rentals",
                type: "bit",
                nullable: true);

            migrationBuilder.AlterColumn<DateTime>(
                name: "PriceDateStart",
                table: "Prices",
                type: "datetime2",
                nullable: false,
                defaultValue: new DateTime(1, 1, 1, 0, 0, 0, 0, DateTimeKind.Unspecified),
                oldClrType: typeof(DateTime),
                oldType: "datetime2",
                oldNullable: true);

            migrationBuilder.AddColumn<double>(
                name: "PricePerNight",
                table: "Prices",
                type: "float",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "RentalId",
                table: "Prices",
                type: "int",
                nullable: false,
                defaultValue: 0);

            migrationBuilder.AlterColumn<int>(
                name: "RentalId",
                table: "Images",
                type: "int",
                nullable: true,
                oldClrType: typeof(int),
                oldType: "int");

            migrationBuilder.AddColumn<bool>(
                name: "PrimaryImage",
                table: "Images",
                type: "bit",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "RentalId",
                table: "Blocked",
                type: "int",
                nullable: true);

            migrationBuilder.CreateTable(
                name: "RentalAmenities",
                columns: table => new
                {
                    Id = table.Column<int>(type: "int", nullable: false)
                        .Annotation("SqlServer:Identity", "1, 1"),
                    RentalId = table.Column<int>(type: "int", nullable: false),
                    AmenityId = table.Column<int>(type: "int", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_RentalAmenities", x => x.Id);
                    table.ForeignKey(
                        name: "FK_RentalAmenities_Amenities_AmenityId",
                        column: x => x.AmenityId,
                        principalTable: "Amenities",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                    table.ForeignKey(
                        name: "FK_RentalAmenities_Rentals_RentalId",
                        column: x => x.RentalId,
                        principalTable: "Rentals",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "RentalBed",
                columns: table => new
                {
                    Id = table.Column<int>(type: "int", nullable: false)
                        .Annotation("SqlServer:Identity", "1, 1"),
                    RentalId = table.Column<int>(type: "int", nullable: false),
                    BedId = table.Column<int>(type: "int", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_RentalBed", x => x.Id);
                    table.ForeignKey(
                        name: "FK_RentalBed_Beds_BedId",
                        column: x => x.BedId,
                        principalTable: "Beds",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                    table.ForeignKey(
                        name: "FK_RentalBed_Rentals_RentalId",
                        column: x => x.RentalId,
                        principalTable: "Rentals",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "Reviews",
                columns: table => new
                {
                    Id = table.Column<int>(type: "int", nullable: false)
                        .Annotation("SqlServer:Identity", "1, 1"),
                    RentalId = table.Column<int>(type: "int", nullable: false),
                    UserName = table.Column<string>(type: "nvarchar(450)", nullable: true),
                    body = table.Column<string>(type: "nvarchar(max)", nullable: false)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_Reviews", x => x.Id);
                    table.ForeignKey(
                        name: "FK_Reviews_AspNetUsers_UserName",
                        column: x => x.UserName,
                        principalTable: "AspNetUsers",
                        principalColumn: "Id");
                    table.ForeignKey(
                        name: "FK_Reviews_Rentals_RentalId",
                        column: x => x.RentalId,
                        principalTable: "Rentals",
                        principalColumn: "Id",
                        onDelete: ReferentialAction.Cascade);
                });

            migrationBuilder.CreateTable(
                name: "UserMails",
                columns: table => new
                {
                    Id = table.Column<int>(type: "int", nullable: false)
                        .Annotation("SqlServer:Identity", "1, 1"),
                    UserName = table.Column<string>(type: "nvarchar(450)", nullable: true),
                    toUser = table.Column<string>(type: "nvarchar(max)", nullable: false),
                    subject = table.Column<string>(type: "nvarchar(max)", nullable: false),
                    body = table.Column<string>(type: "nvarchar(max)", nullable: false),
                    DateSent = table.Column<DateTime>(type: "datetime2", nullable: false),
                    isRead = table.Column<bool>(type: "bit", nullable: true),
                    isDeleted = table.Column<bool>(type: "bit", nullable: true)
                },
                constraints: table =>
                {
                    table.PrimaryKey("PK_UserMails", x => x.Id);
                    table.ForeignKey(
                        name: "FK_UserMails_AspNetUsers_UserName",
                        column: x => x.UserName,
                        principalTable: "AspNetUsers",
                        principalColumn: "Id");
                });

            migrationBuilder.CreateIndex(
                name: "IX_Prices_RentalId",
                table: "Prices",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_Images_RentalId",
                table: "Images",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_Blocked_RentalId",
                table: "Blocked",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_RentalAmenities_AmenityId",
                table: "RentalAmenities",
                column: "AmenityId");

            migrationBuilder.CreateIndex(
                name: "IX_RentalAmenities_RentalId",
                table: "RentalAmenities",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_RentalBed_BedId",
                table: "RentalBed",
                column: "BedId");

            migrationBuilder.CreateIndex(
                name: "IX_RentalBed_RentalId",
                table: "RentalBed",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_Reviews_RentalId",
                table: "Reviews",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_Reviews_UserName",
                table: "Reviews",
                column: "UserName");

            migrationBuilder.CreateIndex(
                name: "IX_UserMails_UserName",
                table: "UserMails",
                column: "UserName");

            migrationBuilder.AddForeignKey(
                name: "FK_Blocked_Rentals_RentalId",
                table: "Blocked",
                column: "RentalId",
                principalTable: "Rentals",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Images_Rentals_RentalId",
                table: "Images",
                column: "RentalId",
                principalTable: "Rentals",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Prices_Rentals_RentalId",
                table: "Prices",
                column: "RentalId",
                principalTable: "Rentals",
                principalColumn: "Id",
                onDelete: ReferentialAction.Cascade);
        }

        /// <inheritdoc />
        protected override void Down(MigrationBuilder migrationBuilder)
        {
            migrationBuilder.DropForeignKey(
                name: "FK_Blocked_Rentals_RentalId",
                table: "Blocked");

            migrationBuilder.DropForeignKey(
                name: "FK_Images_Rentals_RentalId",
                table: "Images");

            migrationBuilder.DropForeignKey(
                name: "FK_Prices_Rentals_RentalId",
                table: "Prices");

            migrationBuilder.DropTable(
                name: "RentalAmenities");

            migrationBuilder.DropTable(
                name: "RentalBed");

            migrationBuilder.DropTable(
                name: "Reviews");

            migrationBuilder.DropTable(
                name: "UserMails");

            migrationBuilder.DropIndex(
                name: "IX_Prices_RentalId",
                table: "Prices");

            migrationBuilder.DropIndex(
                name: "IX_Images_RentalId",
                table: "Images");

            migrationBuilder.DropIndex(
                name: "IX_Blocked_RentalId",
                table: "Blocked");

            migrationBuilder.DropColumn(
                name: "Published",
                table: "Rentals");

            migrationBuilder.DropColumn(
                name: "PricePerNight",
                table: "Prices");

            migrationBuilder.DropColumn(
                name: "RentalId",
                table: "Prices");

            migrationBuilder.DropColumn(
                name: "PrimaryImage",
                table: "Images");

            migrationBuilder.DropColumn(
                name: "RentalId",
                table: "Blocked");

            migrationBuilder.RenameColumn(
                name: "Rental",
                table: "Reservation",
                newName: "RentalId");

            migrationBuilder.RenameColumn(
                name: "Beds",
                table: "Rentals",
                newName: "PriceId");

            migrationBuilder.RenameColumn(
                name: "RentalId",
                table: "Images",
                newName: "RentalID");

            migrationBuilder.AlterColumn<float>(
                name: "SalesTax",
                table: "Reservation",
                type: "real",
                nullable: true,
                oldClrType: typeof(double),
                oldType: "float",
                oldNullable: true);

            migrationBuilder.AlterColumn<float>(
                name: "IndividualTotal",
                table: "Reservation",
                type: "real",
                nullable: true,
                oldClrType: typeof(double),
                oldType: "float",
                oldNullable: true);

            migrationBuilder.AlterColumn<float>(
                name: "GrandTotal",
                table: "Reservation",
                type: "real",
                nullable: false,
                oldClrType: typeof(double),
                oldType: "float");

            migrationBuilder.AlterColumn<string>(
                name: "RentalType",
                table: "Rentals",
                type: "nvarchar(max)",
                nullable: true,
                oldClrType: typeof(string),
                oldType: "nvarchar(max)");

            migrationBuilder.AddColumn<int>(
                name: "AmenityId",
                table: "Rentals",
                type: "int",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "BedId",
                table: "Rentals",
                type: "int",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "BlockId",
                table: "Rentals",
                type: "int",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "ImageId",
                table: "Rentals",
                type: "int",
                nullable: true);

            migrationBuilder.AlterColumn<DateTime>(
                name: "PriceDateStart",
                table: "Prices",
                type: "datetime2",
                nullable: true,
                oldClrType: typeof(DateTime),
                oldType: "datetime2");

            migrationBuilder.AlterColumn<int>(
                name: "RentalID",
                table: "Images",
                type: "int",
                nullable: false,
                defaultValue: 0,
                oldClrType: typeof(int),
                oldType: "int",
                oldNullable: true);

            migrationBuilder.AddColumn<string>(
                name: "Primary",
                table: "Images",
                type: "nvarchar(max)",
                nullable: true);

            migrationBuilder.AddColumn<int>(
                name: "Quantity",
                table: "Beds",
                type: "int",
                nullable: false,
                defaultValue: 0);

            migrationBuilder.CreateIndex(
                name: "IX_Reservation_RentalId",
                table: "Reservation",
                column: "RentalId");

            migrationBuilder.CreateIndex(
                name: "IX_Rentals_AmenityId",
                table: "Rentals",
                column: "AmenityId");

            migrationBuilder.CreateIndex(
                name: "IX_Rentals_BedId",
                table: "Rentals",
                column: "BedId");

            migrationBuilder.CreateIndex(
                name: "IX_Rentals_BlockId",
                table: "Rentals",
                column: "BlockId",
                unique: true,
                filter: "[BlockId] IS NOT NULL");

            migrationBuilder.CreateIndex(
                name: "IX_Rentals_ImageId",
                table: "Rentals",
                column: "ImageId");

            migrationBuilder.CreateIndex(
                name: "IX_Rentals_PriceId",
                table: "Rentals",
                column: "PriceId",
                unique: true,
                filter: "[PriceId] IS NOT NULL");

            migrationBuilder.CreateIndex(
                name: "IX_Images_RentalID",
                table: "Images",
                column: "RentalID",
                unique: true);

            migrationBuilder.AddForeignKey(
                name: "FK_Images_Rentals_RentalID",
                table: "Images",
                column: "RentalID",
                principalTable: "Rentals",
                principalColumn: "Id",
                onDelete: ReferentialAction.Cascade);

            migrationBuilder.AddForeignKey(
                name: "FK_Rentals_Amenities_AmenityId",
                table: "Rentals",
                column: "AmenityId",
                principalTable: "Amenities",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Rentals_Beds_BedId",
                table: "Rentals",
                column: "BedId",
                principalTable: "Beds",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Rentals_Blocked_BlockId",
                table: "Rentals",
                column: "BlockId",
                principalTable: "Blocked",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Rentals_Images_ImageId",
                table: "Rentals",
                column: "ImageId",
                principalTable: "Images",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Rentals_Prices_PriceId",
                table: "Rentals",
                column: "PriceId",
                principalTable: "Prices",
                principalColumn: "Id");

            migrationBuilder.AddForeignKey(
                name: "FK_Reservation_Rentals_RentalId",
                table: "Reservation",
                column: "RentalId",
                principalTable: "Rentals",
                principalColumn: "Id",
                onDelete: ReferentialAction.Cascade);
        }
    }
}
